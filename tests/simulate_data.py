#!/usr/bin/env python3

# Adapted from:
# https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/

import argparse
import collections
import os
import datetime
import numpy as np

parser = argparse.ArgumentParser(
    description='*** Create Read simulation data incorporating amplification' \
        ' and sequencing errors. ***'
)

parser.add_argument(
    '-n', '--number_of_reads', type=int, default=5,
    help='Number of reads. Default = 5.'
)
parser.add_argument(
    '-l', '--length_of_reads', type=int, default=10,
    help='Leangth of the reads sequence (without UMI). Default = 10.'
)
parser.add_argument(
    '-d', '--avg_depth_of_reads', type=int, default=10,
    help='Average read coverage/sequencing depth. Default = 10.'
)
parser.add_argument(
    '-lu', '--length_of_UMIs', type=int, default=5,
    help='Leangth of the UMIs. Default = 5.'
)

# PCR efficiency - guestimate from - http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1298932/
    # initial UMIs, e.g the true number of molecules at a given loci/gene depending on method
parser.add_argument(
    '-emin', '--eff_min', type=float, default=0.8,
    help='Minimum sequencing efficency. Default = 0.8.'
)
parser.add_argument(
    '-emax', '--eff_max', type=float, default=1,
    help='Maximum sequencing efficency. Default = 1.'
)

# Recombinant Taq DNA pol error rate
parser.add_argument(
    '-pcr', '--PCR_cycles', type=int, default=5,
    help='Number of PCR cycles to simulate. Default = 5.'
)
parser.add_argument(
    '-pcr_e', '--PCR_error_rate', type=float, default=0.0001,
    help='PCR amplification error rate. Default = 0.0001.'
)

# Phred = 30
parser.add_argument(
    '-seq_e', '--seq_error_rate', type=float, default=0.001,
    help='Sequencing error rate. Default = 0.001.'
)

parser.add_argument(
    '-o', '--output', type=str, default='',
    help='Absolute or relative path to an output directory.'
)



ERRORS = {
    "A": ["C", "G", "T"],
    "C": ["A", "G", "T"],
    "G": ["C", "A", "T"],
    "T": ["C", "G", "A"]
}


def get_random_seq(n):
    seq = ""
    for base in range(n):
        seq += np.random.choice(["A", "C", "G", "T"])
    return seq


def createReadsWithBias(n, read_size, read_depth, UMI_size, efficiency_min,
            efficiency_max):

    assert 0 < efficiency_min, "PCR efficiency minimum must be between 0 and 1"
    assert efficiency_min <= 1, "PCR efficiency maximum must be between 0 and 1"

    reads = []
    read_to_UMI = {}
    efficiency_dict = collections.defaultdict(float)
    depth_dict = {}
    for read in range(n):
        seq = get_random_seq(read_size)
        depth = np.random.poisson(read_depth)
        depth_dict[seq] = depth

        for i in range(depth):
            UMI = get_random_seq(UMI_size)
            read = UMI + seq
            reads.append(read)
            efficiency_dict[read] = np.random.uniform(
                efficiency_min, efficiency_max
            )
            try:
                read_to_UMI[seq].append(UMI)
            except KeyError:
                read_to_UMI[seq] = [UMI]

    return reads, efficiency_dict, read_to_UMI, depth_dict


def PCRcyclesWithErrorsAndBias(reads, efficiency, PCR_cycles, error_rate):
    for cycle in range(PCR_cycles):
        post_cycle = simulatePCRcycleWithErrorsAndBias(
            reads, efficiency, error_rate
        )
        reads = post_cycle
    return post_cycle


def simulatePCRcycleWithErrorsAndBias(reads, efficiency, error_rate):
    new_list = []
    for read in reads:
        new_list.extend(
            amplifyWithErrors(read, efficiency[read], error_rate)
        )
    return new_list


def amplifyWithErrors(read, efficiency, error_rate):
    '''simulate amplification on a single UMI
    efficiency'''
    if np.random.random() > (1 - efficiency):
        if error_rate == 0:
            return (read, read)
        else:
            err_read = get_err_read(read, error_rate)
            return (read, err_read)
    else:
        return (read,)


def get_err_read(read, error):
    err_read = ""
    for base in read:
        if np.random.random() <= error:
            err_read += np.random.choice(ERRORS[base])
        else:
            err_read += base
    return err_read


def addSequencingErrors(reads, seq_error_rate):
    '''takes a UMI counter object and adds errors at random to simulate sequencing errors'''
    err_reads = [get_err_read(i, seq_error_rate) for i in reads]
    return err_reads


def save(obj, out_dir, file_name):
    out_file = os.path.join(out_dir, file_name)

    if isinstance(obj, list):
        _save_list(obj, out_file)
    elif isinstance(obj, dict):
        _save_dict(obj, out_file)
    else:
        raise IOError('Unsupported save type: {}'.format(type(obj)))


def _save_list(reads, file):
    with open(file, 'w') as f:
        f.write('\n'.join(reads))


def _save_dict(out_dict, file):
    with open(file, 'w') as f:
        for key, val in out_dict.items():
            f.write('{}\t{}\n'.format(key, val))


def save_as_fastq(reads, out_dir, file_name):
    out_file = os.path.join(out_dir, file_name)
    fastq_str = ''
    for idx, read in enumerate(reads):
        read_str = '@SEQ_{}\n{}\n+\n{}\n'.format(idx, read, '*' * len(read))
        fastq_str += read_str
    with open(out_file, 'w') as f:
        f.write(fastq_str)


if __name__ == '__main__':
    args = parser.parse_args()

    reads, efficiency, read_to_UMI, depth = createReadsWithBias(
        args.number_of_reads, args.length_of_reads, args.avg_depth_of_reads,
        args.length_of_UMIs, args.eff_min, args.eff_max
    )

    PCR_reads = PCRcyclesWithErrorsAndBias(
        reads, efficiency, args.PCR_cycles, args.PCR_error_rate
    )

    seq_reads = addSequencingErrors(PCR_reads, args.seq_error_rate)

    if args.output:
        out_dir = args.output
    else:
        out_dir = 'simulations'

    if args.output:
        out_dir = args.output
    else:
        out_dir = '{:%Y%m%d_%H%M%S}'.format(datetime.datetime.now())

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    save(vars(args), out_dir, 'args.txt')

    save(read_to_UMI, out_dir, 'seq_UMI_map.txt')
    save(reads, out_dir, 'reads_raw.txt')
    save_as_fastq(seq_reads, out_dir, 'reads_raw.fq')
    save(collections.Counter(reads), out_dir, 'reads_raw_counts.txt')
    save(depth, out_dir, 'results.txt')

    save(PCR_reads, out_dir, 'reads_PCR.txt')
    save(collections.Counter(PCR_reads), out_dir, 'reads_PCR_counts.txt')
    save_as_fastq(seq_reads, out_dir, 'reads_PCR_counts.fq')


    save(seq_reads, out_dir, 'reads_PCR_Seq.txt')
    save(collections.Counter(seq_reads), out_dir, 'reads_seq_counts.txt')
    save_as_fastq(seq_reads, out_dir, 'reads_PCR_Seq.fq')

    print("number of reads before sequencing: {}".format(len(reads)))
    print("number of reads after PCR: {}".format(len(PCR_reads)))
    print("number of reads after sequencing: {}".format(len(seq_reads)))

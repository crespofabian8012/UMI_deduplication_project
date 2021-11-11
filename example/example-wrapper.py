#!/usr/bin/env python3

import argparse
import os
from datetime import datetime

from umidedup import UMIset


# ------------------------------------------------------------------------------
# ARGPARSER
# ------------------------------------------------------------------------------

def check_fastq(file_path):
    if not file_path.endswith(('fq', '.fastq')):
        raise TypeError(
            'Invalid fastqq file format. Supported formats are: fq|fastq'
        )
    if not os.path.isfile(file_path):
        raise FileNotFoundError()
    return file_path


def check_output(output_path):
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    return output_path


parser = argparse.ArgumentParser(
    description='*** PCR and Sequencing Error correction for Unique molecular '
        'identifiers (UMIs) with subsequent deduplication. ***'
)

parser.add_argument(
    '-i', '--input', type=check_fastq, required=True,
    help='Absolute or relative path to the input data in fastq file format.'
)

parser.add_argument(
    '-l', '--length', type=int, default=7,
    help='Length of UMIs in base pairs. Default = 7.'
)

parser.add_argument(
    '-a', '--algorithm', type=str,
    choices=['hamming', 'alevin', 'frequency'],
    default='alevin', help='Algorithm used for UMI deduplication. ' \
        'Default = "alevin"\n\n' \
    '- Hamming: Compute pairwise hamming distance between UMI sequences and merge those UMIs that have a hamming' \
    'distance smaller or equal to the threshold value, the merging is done in a greedy manner, starting and growing' \
    'those UMIs that have the most reads, no reads are discarded'
    '- Frequency: Remove UMIs with less than threshold reads from the fastq file'
)

parser.add_argument(
    '-t', '--threshold', type=int, default=-1,
    help='Cutoff treshold for algorithms "hamming" and "frequency". ' \
        'Default for Hamming = 1, Default for frequency = 2' \
)

parser.add_argument(
    '-o', '--output', type=check_output, default='UMI_output',
    help='Absolute or relative path to an output directory. The processed fastq' \
        ' files as well as UMI error correction and UMI deduplication statistics '\
        'will be saved here. Default = "./UMI_output"'
)

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------

if __name__ == '__main__':
    start_time = datetime.now()
    print('\nStart time UMI deduplication:\t{}'.format(start_time))

    args = parser.parse_args()
    print('\tInput file:\t{}'.format(args.input))
    print('\tAlgorithm:\t{}'.format(args.algorithm))
    print('\tOutput dir:\t{}\n'.format(args.output))

    print('1. Reading UMIs')
    umi_set = UMIset.UMIset(args.input, args.length)
    print('1. Reading UMIs - done')

    print('2. Deduplicating UMIs')
    if args.threshold == -1:
        if args.algorithm == 'hamming':
            args.threshold = 1
        elif args.algorithm == 'frequency':
            args.threshold = 2
        elif args.algorithm == 'alevin':
            args.threshold = 3
    umi_set.collapse(method=args.algorithm, threshold=args.threshold)
    print('2. Deduplicating UMIs - done')

    print('3. Saving results')
    umi_set.save_statistics(args.output)
    umi_set.save_processed_fastq(args.output)
    print('3. Saving results - done')

    end_time = datetime.now()
    print('End time UMI collapsing:\t{}'.format(end_time))
    print('\tTotal time elapsed: {}'.format(end_time - start_time))

#!/usr/bin/env python3

import os
import sys
import tempfile
import pytest
from umidedup import fastqIO
from umidedup import UMIset

BASE_DIR = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-2])
sys.path.insert(0, os.path.join(BASE_DIR, 'umidedup'))


fastq = os.path.join(BASE_DIR, 'tests', 'test1.fq')


data_write_fastq = ([fastq])
@pytest.mark.parametrize('fastq_file', data_write_fastq)
def test_write_fastq(fastq_file):
    data = fastqIO.read(fastq_file)
    out_file = tempfile.NamedTemporaryFile(delete=False)
    fastqIO.write(data, out_file.name)

    data2 = fastqIO.read(out_file.name)
    assert sorted(data.keys()) == sorted(data2.keys())
    assert sorted(data.values()) == sorted(data2.values())


data_read_fastq = [
    (fastq, 'AGGGGC', 0, 'GTTCGTTCAAGTGCACTTTCCAGTACACTTA'),
    (fastq, 'CAGAGC', 0, 'CCCTTTTTCCCCCAGATCGGAAAAACACACC'),
    (fastq, 'TTAAGG', 0, 'GCTACTACCACCAAGATCTGCACCTGCGGCG'),
    (fastq, 'TTAAGG', 1, 'AGGGTGGGGGATCACATTTATTGTATTGAGG')
]
@pytest.mark.parametrize('fastq_file, UMI, no, expected', data_read_fastq)
def test_read_fastq(fastq_file, UMI, no, expected):
    data = fastqIO.read(fastq_file)
    assert data[UMI][no].get_seq() == expected


@pytest.mark.parametrize('fastq_file, UMI, no, expected', data_read_fastq)
def test_write_Umiset(fastq_file, UMI, no, expected):
    # Init and collapse
    umiset = UMIset.UMIset(fastq_file, 6)
    umiset.collapse('frequency')
    # Save output
    out_file = tempfile.NamedTemporaryFile(delete=False)
    umiset.save_processed_fastq(out_file.name)
    assert fastqIO.read(fastq_file) == fastqIO.read(out_file.name)

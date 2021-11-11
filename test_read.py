#!/usr/bin/env python3

import os
import sys
import pytest
from umidedup.read import Read
from umidedup.UMIset import UMIset

BASE_DIR = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-2])
sys.path.insert(0, os.path.join(BASE_DIR, 'umidedup'))

fastq = os.path.join(BASE_DIR, 'tests', 'test1.fq')

fastq_data = [
    ('@id1', 'AAAAACCCCC', '*****'), ('@id2', 'aaaaacccccc', '*****'),
    ('@id3', 'ACGTACGTA', '*****'), ('@id4', 'GGGGGGGGGGGGGGGGGGGGGG', '*****'),
]

@pytest.mark.parametrize('id_str, seq, qual', fastq_data)
def test_id(id_str, seq, qual):
    read = Read(id_str, seq, qual)
    assert read.get_id() == id_str


@pytest.mark.parametrize('id_str, seq, qual', fastq_data)
def test_seq(id_str, seq, qual):
    read = Read(id_str, seq, qual)
    assert read.get_seq() == seq.upper()


@pytest.mark.parametrize('id_str, seq, qual', fastq_data)
def test_cell_ID_trimming(id_str, seq, qual):
    read = Read(id_str, seq, qual)
    bp = 3
    read.trim_cell_ID(bp)
    assert read.get_cell_ID() == seq[:bp].upper()


@pytest.mark.parametrize('id_str, seq, qual', fastq_data)
def test_UMI_trimming(id_str, seq, qual):
    read = Read(id_str, seq, qual)
    bp = 5
    read.trim_UMI(bp)
    assert read.get_UMI() == seq[:bp].upper()

#!/usr/bin/env python3

import os
import sys
import pytest
from umidedup import fastqIO
from umidedup import UMIset
from umidedup import collapsing

BASE_DIR = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-2])
sys.path.insert(0, os.path.join(BASE_DIR, 'src'))

@pytest.fixture
def fastq_fixture():
    fastq_file = os.path.join(BASE_DIR, 'tests', 'test_collapse.fq')
    return fastq_file

def test_collapse_hamming(fastq_fixture):
    umiset = UMIset.UMIset(fastq_fixture, 6)
    umap = collapsing.collapse_hamming(umiset, threshold=3)
    assert(list(umap.keys()) == ["AAAAAA", "TTTTTT"])

def test_collapse_frequency(fastq_fixture):
    umiset = UMIset.UMIset(fastq_fixture, 6)
    umap = collapsing.collapse_frequency(umiset, min_count=3)
    test = []
    for k,v in umap.items():
        if v == []:
            test.append(k)

    assert(test == ["AAAAAC", "TTTTTA", "AAAAAG"])

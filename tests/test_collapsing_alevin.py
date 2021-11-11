
import os
import sys
import pytest
import importlib

BASE_DIR = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-2])
sys.path.insert(0, os.path.join(BASE_DIR, 'umidedup'))

from collections import defaultdict, OrderedDict
from itertools import combinations

from umidedup.UMIset import UMIset
from umidedup.collapsing import build_network
from umidedup.collapsing import get_all_connected_components
from umidedup.collapsing import process_connected_components

fastq = os.path.join(BASE_DIR, 'tests', 'test3.fq')

data_read_fastq = [
    (fastq, 'AGGGGC', 0, 'GTTCGTTCAAGTGCACTTTCCAGTACACTTA'),
    (fastq, 'CAGAGC', 0, 'CCCTTTTTCCCCCAGATCGGAAAAACACACC'),
    (fastq, 'TTAAGG', 0, 'GCTACTACCACCAAGATCTGCACCTGCGGCG'),
    (fastq, 'TTAAGG', 1, 'AGGGTGGGGGATCACATTTATTGTATTGAGG')
]
@pytest.mark.parametrize('fastq_file, UMI, no, expected', data_read_fastq)
def test_build_network(fastq_file, UMI, no, expected):
    Umiset= UMIset(fastq_file, 6)
    UMis= Umiset.get_UMIs()
    UMInetwork = build_network(Umiset)
    assert UMI  in UMInetwork.keys()
    assert sorted(UMInetwork['ATTCCG']) == sorted({'ATACCG', 'AATCCG', 'ATTCCA'})

fastq = os.path.join(BASE_DIR, 'tests', 'test3.fq')

data_read_fastq = [
    (fastq, 'AGGGGC', 0, 'GTTCGTTCAAGTGCACTTTCCAGTACACTTA'),
    (fastq, 'CAGAGC', 0, 'CCCTTTTTCCCCCAGATCGGAAAAACACACC'),
    (fastq, 'TTAAGG', 0, 'GCTACTACCACCAAGATCTGCACCTGCGGCG'),
    (fastq, 'TTAAGG', 1, 'AGGGTGGGGGATCACATTTATTGTATTGAGG')
]
@pytest.mark.parametrize('fastq_file, UMI, no, expected', data_read_fastq)
def test_get_all_connected_components(fastq_file, UMI, no, expected):
    Umiset = UMIset(fastq_file, 6)
    UMis = Umiset.get_UMIs()
    UMInetwork = build_network(Umiset)
    list_connected_components = get_all_connected_components(UMInetwork)
    assert len(list_connected_components)==230

fastq = os.path.join(BASE_DIR, 'tests', 'test3.fq')

data_read_fastq = [
    (fastq, 'AGGGGC', 0, 'GTTCGTTCAAGTGCACTTTCCAGTACACTTA'),
    (fastq, 'CAGAGC', 0, 'CCCTTTTTCCCCCAGATCGGAAAAACACACC'),
    (fastq, 'TTAAGG', 0, 'GCTACTACCACCAAGATCTGCACCTGCGGCG'),
    (fastq, 'TTAAGG', 1, 'AGGGTGGGGGATCACATTTATTGTATTGAGG')
]
@pytest.mark.parametrize('fastq_file, UMI, no, expected', data_read_fastq)
def test_process_connected_components(fastq_file, UMI, no, expected):
    Umiset = UMIset(fastq_file, 6)
    UMis = Umiset.get_UMIs()
    UMInetwork = build_network(Umiset)
    list_connected_components = get_all_connected_components(UMInetwork)
    dict_collapsing_UMIs = process_connected_components(UMInetwork, Umiset)
    assert len(dict_collapsing_UMIs)== 230

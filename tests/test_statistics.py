
import os
import sys
from collections import defaultdict, OrderedDict
from itertools import combinations
import pytest
import importlib

from umidedup.UMIset import UMIset
from umidedup.collapsing import build_network
from umidedup.collapsing import get_all_connected_components
from umidedup.collapsing import process_connected_components
from umidedup import statistics

BASE_DIR = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-2])
sys.path.insert(0, os.path.join(BASE_DIR, 'umidedup'))

fastq = os.path.join(BASE_DIR, 'example','example_data', 'reads.fq')
ground_truth = os.path.join(BASE_DIR, 'example_data', 'ground_truth.txt')

ground_truth_dict = {
'TATAGTATGAAGGTAGTGAAGCTGACTTCCGACCCCCTTACTTTTCCCAG':97,
'CTACGGTAGGAAAACTATGGTGGCACTCCTGTTATGAAACTGACGCCAAC':93,
'GGAATGTATAGCCTCCATCAAAAGGGCCGTATGTTACGAGTCTAATCGAT':111,
'TGATTAAATTGAGAAACGTATCGTTTCCCAGATTACGGGTAAGCGTATAG':100,
'CCCGTCTAGTTCTCTCAAACACCGCATGTAGGATTGTCGTGGGACCAATG':73,
'CCCGGGTATCGATCGGATCTCTCTTTCTCTTCGCCCCGGTTACTCCGTCT':111,
'CTGACCTCTTTACCTTACCACCTTAGGGGGAATTAGTGATGGACCAGCAG':122,
'ATGTGATATTCGGTTGCACCGCCTGTCCTAGCCATATTTGGCAACCTTAG':110,
'CACGACCGTGTCCCTACACGACAAGTTAAACCGCTTTCCTAGTGCGGGCC':103,
'TTCGCGTGGTGCGTGCAGACACTGTAGGAGTCCACAATAATTATGACTAG':80
}

data_compute_statistics= [
    (fastq)
]
@pytest.mark.parametrize('fastq_file', data_compute_statistics)
def test_compute_statistics_one_simulation(fastq_file):
    Umiset = UMIset(fastq_file, 7)
    Umiset.collapse( "alevin",  1)
    (number_collapsed_UMIs, relative_number_collapsed_UMIs, relative_number_unique_UMIs, \
     number_unique_UMIs_before_collapsing, \
     avg_number_of_reads_before_collapsing, number_of_collapsed_reads, proportion_of_collapsed_reads,number_of_reads_before_collapsing,avg_number_of_reads_after_collapsing)=   statistics.compute_statistics_one_simulation(Umiset)
    assert number_of_collapsed_reads == 839
    assert number_collapsed_UMIs == 152
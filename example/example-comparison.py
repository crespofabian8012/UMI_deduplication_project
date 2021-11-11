#!/usr/bin/env python3

import os
import sys
from umidedup.UMIset import UMIset
from umidedup import statistics

BASE_DIR = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-2])

def comparison(fastqpath, umilength):

    print("Applying naive frequency collapsing")
    umiset1 = UMIset(filepath, 7)
    umiset1.collapse("frequency", threshold = 10)
    (number_collapsed_UMIs, relative_number_collapsed_UMIs, relative_number_unique_UMIs, \
     number_unique_UMIs_before_collapsing, \
     avg_number_of_reads_before_collapsing, number_of_collapsed_reads, proportion_of_collapsed_reads,number_of_reads_before_collapsing,avg_number_of_reads_after_collapsing) = statistics.compute_statistics_one_simulation(umiset1)

    print("Statistics before collapsing")
    print("# UMIs: ",number_unique_UMIs_before_collapsing )
    print("# reads: ", number_of_reads_before_collapsing)
    print("Average #reads per UMI: ", avg_number_of_reads_before_collapsing)

    print("Statistics after collapsing")
    print("# collapsed UMIs: ", number_collapsed_UMIs, " relative # of collapsed UMIs:", relative_number_collapsed_UMIs)
    print("# collapsed reads: ", number_of_collapsed_reads,
          "Proportion of #collapsed reads to the total number of reads: ", proportion_of_collapsed_reads)
    print("Average #reads per UMI: ", avg_number_of_reads_after_collapsing)


    print()
    print("Applying hamming collapsing")
    umiset2 = UMIset(filepath, 7)
    print(len(umiset2.UMIs.keys()))
    umiset2.collapse("hamming", 3)

    (number_collapsed_UMIs2, relative_number_collapsed_UMIs2, relative_number_unique_UMIs2, \
     number_unique_UMIs_before_collapsing2, \
     avg_number_of_reads_before_collapsing2, number_of_collapsed_reads2,\
     proportion_of_collapsed_reads2,number_of_reads_before_collapsing2,avg_number_of_reads_after_collapsing2) = statistics.compute_statistics_one_simulation(umiset2)


    print("Statistics before collapsing")
    print("# UMIs: ", number_unique_UMIs_before_collapsing2)
    print("# reads: ", number_of_reads_before_collapsing2)

    print("Average #reads per UMI: ", avg_number_of_reads_before_collapsing2)

    print("Statistics after collapsing")
    print("# collapsed UMIs: ", number_collapsed_UMIs2, " relative # of collapsed UMIs:", relative_number_collapsed_UMIs2)
    print("# collapsed reads: ", number_of_collapsed_reads2,"Proportion of #collapsed reads to the total number of reads: ", proportion_of_collapsed_reads2)
    print("Average #reads per UMI: ", avg_number_of_reads_after_collapsing2)

    print()
    print("Applying Alevin collapsing")
    umiset3 = UMIset(filepath, 7)
    umiset3.collapse("alevin", 3)

    (number_collapsed_UMIs3, relative_number_collapsed_UMIs3, relative_number_unique_UMIs3, \
     number_unique_UMIs_before_collapsing3, \
     avg_number_of_reads_before_collapsing3, number_of_collapsed_reads3,\
     proportion_of_collapsed_reads3,number_of_reads_before_collapsing3,avg_number_of_reads_after_collapsing3) = statistics.compute_statistics_one_simulation(umiset3)

    print("Statistics before collapsing")
    print("# UMIs: ", number_unique_UMIs_before_collapsing3)
    print("# reads: ", number_of_reads_before_collapsing3)
    print("Average #reads per UMI: ", avg_number_of_reads_before_collapsing3)

    print("Statistics after collapsing")
    print("# collapsed UMIs: ", number_collapsed_UMIs3, " relative # of collapsed UMIs:", relative_number_collapsed_UMIs3)
    print("# collapsed reads: ", number_of_collapsed_reads3,
          "Proportion of #collapsed reads to the total number of reads: ", proportion_of_collapsed_reads3)
    print("Average #reads per UMI: ", avg_number_of_reads_after_collapsing3)

if __name__ == "__main__":
    print("Example using reads.fq file in example_data folder")
    filepath = os.path.join(BASE_DIR, "example", "example_data", "reads.fq")
    comparison(filepath, 7)

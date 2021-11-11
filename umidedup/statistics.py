#!/usr/bin/env python3

from umidedup import UMIset
from umidedup import UMImapping
import statistics as stat

def compute_statistics_one_simulation(UMIset):
    if not UMIset.UMIs_processed:
        raise RuntimeError('No processed UMIs found')
    #umap = UMIset.umap
    #dict_number_collapsed_UMIs= {key: len(value) for key, value in umap.items()}
   # number_collapsed_UMIs = sum(dict_number_collapsed_UMIs.values())
    number_collapsed_UMIs = len(UMIset.UMIs )-len(UMIset.UMIs_processed )
    relative_number_collapsed_UMIs = number_collapsed_UMIs \
        / len(UMIset.get_list_UMIs())
    number_unique_UMIs_before_collapsing= len(UMIset.get_list_UMIs())

    relative_number_unique_UMIs = number_unique_UMIs_before_collapsing \
        / len(UMIset.get_list_UMIs())
    dict_number_reads_before_collapsing = {
        UMI: len(UMIset.get_list_reads_for_UMI(UMI)) \
            for  UMI in set(UMIset.get_list_UMIs())
    }
    avg_number_of_reads_before_collapsing = sum(
        dict_number_reads_before_collapsing.values()) \
            / len(UMIset.get_list_UMIs())

    #umimapping = {z: x for x, y in umap.items() for z in y}
    # number_of_collapsed_reads= sum(dict_number_reads_before_collapsing[UMI] for UMI in set(umimapping.keys()))
    collapsed_UMis = set(set(UMIset.get_list_UMIs()) - set(UMIset.UMIs_processed.keys()))
    number_of_collapsed_reads =  sum(dict_number_reads_before_collapsing[UMI] for UMI in collapsed_UMis)
    number_of_reads_before_collapsing = sum(
        dict_number_reads_before_collapsing[UMI] for UMI in set(UMIset.get_list_UMIs()))

    proportion_of_collapsed_reads = number_of_collapsed_reads \
        / sum(dict_number_reads_before_collapsing.values())

    avg_number_of_reads_after_collapsing =  (number_of_reads_before_collapsing - number_of_collapsed_reads) / ( number_unique_UMIs_before_collapsing - number_collapsed_UMIs)

    return (number_collapsed_UMIs, relative_number_collapsed_UMIs, relative_number_unique_UMIs, number_unique_UMIs_before_collapsing, \
        avg_number_of_reads_before_collapsing, number_of_collapsed_reads, proportion_of_collapsed_reads, number_of_reads_before_collapsing,avg_number_of_reads_after_collapsing)

def compute_statistics_for_multiple_simulations(list_UMI_mappings, list_UMISet):
    list_relative_number_unique_UMIs =[]
    list_relative_number_collapsed_UMIs  =[]

    for pair in zip(list_UMI_mappings, list_UMISet):
        UmiMapping=pair[0]
        UmiSet=pair[1]
        relative_number_unique_UMIs, relative_number_collapsed_UMIs = \
        compute_statistics_one_simulation(UmiMapping, UmiSet)
        list_relative_number_unique_UMIs.append(relative_number_unique_UMIs)
        list_relative_number_collapsed_UMIs.append(relative_number_collapsed_UMIs)

    coef_var_relative_number_unique_UMIs = \
    stat.stdev(list_relative_number_unique_UMIs) \
    / stat.mean(list_relative_number_unique_UMIs)
    coef_variation_relative_number_collapsed_UMIs = \
    stat.stdev(list_relative_number_collapsed_UMIs) \
    / stat.mean(list_relative_number_collapsed_UMIs)
    return (coef_var_relative_number_unique_UMIs,
        coef_variation_relative_number_collapsed_UMIs)

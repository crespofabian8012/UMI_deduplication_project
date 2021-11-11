#!/usr/bin/env python3

from collections import defaultdict, OrderedDict
from itertools import combinations
from umidedup import UMImapping
from umidedup import UMIset

def hamming_distance(UMI0, UMI1):
    """Hamming distance between UMIs (assumes equal length)"""
    if len(UMI0) != len(UMI1):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(l1 != l2 for l1, l2 in zip(UMI0, UMI1))

def collapse_frequency(UMIset:UMIset, min_count=3):
    """ remove UMIs with less or equal to min_count reads """
    umap = UMImapping.UMImapping().umap
    UMIs =  UMIset.get_UMIs()
    for (key, value) in UMIs:
        if len(value) <= min_count:
            umap[key] = []
        else:
            umap[key] = [key]

    return umap

def collapse_hamming(UMIset:UMIset, threshold=3):
    """ simple merging based on hamming distance,
        greedy process start with UMIs with highest count to avoid ambiguity """

    # get counts for each UID
    counts = OrderedDict()
    UMIs = UMIset.get_UMIs()
    umap = UMImapping.UMImapping().umap

    UMI_seqs = []
    for (key, value) in UMIs:

        umap[key] = [key]
        if isinstance(value, int):
            val = value
        else:
            val = len(value)

        if val in counts:
            counts[val].append(key)
        else:
            counts[val] = [key]

    # iterate over sorted counts and merge greedily, starting with largest
    sorted_counts = sorted(counts)[::-1]

    for c in sorted_counts:
        try:
            umi_cs = counts[c]
        except KeyError:
            raise("Key error")

        for umi_c in umi_cs:

            tmp_keys = list(umap.keys())
            if umi_c not in tmp_keys:
                continue

            for umi_other in tmp_keys:

                if umi_c == umi_other:
                    continue

                dist = hamming_distance(umi_c, umi_other)

                if dist <= threshold:
                    # update mapping
                    if umi_c in umap:
                        umap[umi_c].extend(umap[umi_other])
                    if umi_other in umap:
                        del umap[umi_other]

    return umap

### Alevin methods ###
def build_network(UMIset):
        UMInetwork = defaultdict(int)
        allpairsUMIs = list(combinations(UMIset.get_list_UMIs(), 2))
        distance = 0
        for pairUMI in allpairsUMIs:
            if pairUMI[0] not in UMInetwork.keys():
               UMInetwork[pairUMI[0]]=set()
            if pairUMI[1] not in UMInetwork.keys():
                UMInetwork[pairUMI[1]] = set()
            distance = hamming_distance(pairUMI[0], pairUMI[1])
            len_intersection=len(intersection_set_reads(pairUMI[0], pairUMI[1], UMIset))
            if distance == 0 and len_intersection > 0:
                _add_edge(pairUMI[0], pairUMI[1], UMInetwork, UMIset)
                _add_edge(pairUMI[1], pairUMI[0], UMInetwork, UMIset)
            elif distance == 1 and len_intersection > 0:
                if _possible_amplification_error(pairUMI[0],  pairUMI[1], UMIset):
                    _add_edge(pairUMI[0], pairUMI[1], UMInetwork, UMIset)
                elif _possible_amplification_error( pairUMI[1],  pairUMI[0], UMIset):
                    _add_edge( pairUMI[1], pairUMI[0], UMInetwork, UMIset)

        return UMInetwork

def _add_edge(UMI1, UMI2, UMInetwork, UMIset):
    if UMI1 in UMIset.get_list_UMIs():
        UMInetwork[UMI1].add(UMI2)
    else:
        UMInetwork[UMI1] = set([UMI2])

def intersection_set_reads(UMI1, UMI2,  UMIset):
    UMIs=UMIset.get_list_UMIs()
    if UMI1 in UMIs and UMI2 in UMIs:
        s1 = set(UMIset.get_list_reads_for_UMI(UMI1))#this is a set
        s2 = set(UMIset.get_list_reads_for_UMI(UMI2))#this is another set
        return s1.intersection(s2)
    else:
       return set()

def  cardinality_set(set):
       return len(set)

def _add_bireccional_edge(UMI1, UMI2, UMInetwork, UMIset):
    _add_edge(UMI1, UMI2, UMInetwork, UMIset)
    _add_edge(UMI2, UMI1, UMInetwork, UMIset)

def _possible_pcr_error(UMI1, UMI2, UMIset):
    return len(set(UMIset.get_list_reads_for_UMI(UMI1))) > 2 * len(set(UMIset.get_list_reads_for_UMI(UMI2))) - 1

def _possible_amplification_error(UMI1, UMI2, UMIset):
    return _possible_pcr_error(UMI1, UMI2, UMIset)

def depth_first_search_from_vertex(UMInetwork, start, already_visited):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(UMInetwork[vertex] - visited)
    return visited

def get_all_connected_components(UMInetwork):
    already_visited = set()
    result = []
    for node in UMInetwork.keys():
        if node not in already_visited:
            connected_group, already_visited = get_connected_component_from_node(UMInetwork, node, already_visited)
            result.append(connected_group)
    return result

def get_connected_component_from_node(UMInetwork, node, already_visited):
        result = []
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            already_visited.add(node)
            nodes = nodes or UMInetwork[node] - already_visited
            result.append(node)
        return result, already_visited

def breath_first_search_from_node(UMInetwork, start):
    visited, queue = set(), [start]
    while queue:
        vertex = queue.pop(0)
        if vertex not in visited:
            visited.add(vertex)
            queue.extend(UMInetwork[vertex] - visited)
    return visited

def process_connected_components(UMInetwork, UMIset):
    list_connected_components= get_all_connected_components(UMInetwork)
    dict_collapsing_UMIs = OrderedDict()
    for component in list_connected_components:
        length_tree_from_node = {start: len(breath_first_search_from_node(UMInetwork, start)) for start  in component}
        UMI_largest_tree = max(length_tree_from_node.keys(), key=(lambda k: length_tree_from_node[k]))
        dict_collapsing_UMIs[UMI_largest_tree]= set(component) - set([UMI_largest_tree])
    return dict_collapsing_UMIs

def collapse_alevin(UMIset:UMIset):
    umap = OrderedDict()
    UMIs = UMIset.get_UMIs()
    try:
       UMInetwork = build_network(UMIset)
       umap = process_connected_components(UMInetwork, UMIset)
    except Exception as e:
        raise

    return umap

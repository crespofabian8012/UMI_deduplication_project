#!/usr/bin/env python3

import os
import re
from collections import OrderedDict
from umidedup.read import Read

def read(in_file, UMI_bp=5, cell_ID_bp=0):
    """Reads a fastq file into a dictionary where keys are UMIs and values are
        reads

    Args:
        in_file (str): Absolute or relative path to FASTq file
        UMI_bp (int): Number of bais pairs encoding each UMI. Default = 5.
        cell_ID_bp (int): Number of bais pairs encoding each cell ID. Default = 0.

    Returns:
        dict (str -> list of Read): Keys are UMIs, values are reads with the
            corresponding UMI
    """
    #  data_name_raw = os.path.basename(in_file)
    with open(in_file, 'r') as f:
        data_raw = f.read().strip().split('\n')

    UMIs = OrderedDict()
    for i in range(0, len(data_raw), 4):

        new_read = Read(
            ident=data_raw[i],
            seq=data_raw[i+1],
            qual=data_raw[i+3]
        )
        if cell_ID_bp:
            new_read.trim_cell_ID(cell_ID_bp)

        if 'UMI_' in data_raw[i]:
            new_UMI = re.search('UMI_([A|C|G|T]+)', data_raw[i]).group(1)
            new_read.set_UMI(new_UMI)
        else:
            new_UMI = new_read.trim_UMI(UMI_bp)

        try:
            UMIs[new_UMI].append(new_read)
        except KeyError:
            UMIs[new_UMI] = [new_read]

    return UMIs


def write(UMIs, out_file):
    """Writes a set of deduplicated reads to a fastq file

    Args:
        UMIs dict|list: Deduplicated reads either as a list or as dict values.

    Returns:
        None.
    """
    if isinstance(UMIs, dict):
        UMIs = [i for j in UMIs.values() for i in j]

    with open(out_file, 'w') as f:
        for r in UMIs:
            f.write(str(r))

#!/usr/bin/env python3

from umidedup import fastqIO
from umidedup import collapsing

# ------------------------------------------------------------------------------
# UMIset
# ------------------------------------------------------------------------------


class UMIset:
    def __init__(self, fastq_file, UMI_length):
        self.UMIs = self._read_fastq(fastq_file, UMI_length)
        self.UMIs_processed = None
        self.umap = None

    def get_UMIs(self):
        return self.UMIs.items()

    def get_UMIs_dict(self):
        return self.UMIs

    def get_list_reads_for_UMI(self, UMI):
        if UMI in self.UMIs:
            return [read.get_seq() for read in self.UMIs[UMI]]

    def get_list_UMIs(self):
        return self.UMIs.keys()

    def _read_fastq(self, infile, UMI_length):
        return fastqIO.read(infile, UMI_bp=UMI_length)

    def save_processed_fastq(self, outfile):
        if self.UMIs_processed:
            return fastqIO.write(self.UMIs_processed, outfile)
        else:
            raise RuntimeError('No processed UMIs found')

    def collapse(self, method, threshold=3):

        # call method and apply to UMIset
        if method == "frequency":
           self.umap = collapsing.collapse_frequency(self, min_count=threshold)
        elif method == "hamming":
            self.umap = collapsing.collapse_hamming(self, threshold=threshold)
        elif method == "alevin":
            self.umap = collapsing.collapse_alevin(self)
        else:
            raise("Method not supported")

        # apply umap to UMIset
        umimapping = {z: x for x, y in self.umap.items() for z in y}

        new_UMIs = dict()
        for k, v in self.UMIs.items():
            try:
                new_UMIs[umimapping[k]] = self.UMIs[k]
            except KeyError:
                if not self.umap[k] == []:
                    new_UMIs[k] = v

        self.UMIs_processed = new_UMIs

# ------------------------------------------------------------------------------

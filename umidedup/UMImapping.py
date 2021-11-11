#!/usr/bin/env python3

from collections import OrderedDict
# ------------------------------------------------------------------------------
# UMImapping
# ------------------------------------------------------------------------------

class UMImapping:

    def __init__(self):
        # object is a dictionary of the form
        # Dict(UMI -> list of UMI)
        self.umap= OrderedDict()

# ------------------------------------------------------------------------------

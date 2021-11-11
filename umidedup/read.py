#!/usr/bin/env python3

# ------------------------------------------------------------------------------
# Read
# ------------------------------------------------------------------------------


class Read:
    """ Representation of a sequenced read in FASTQ format.

    A read consists of three elements:
        1. A identifier, starting with an "@" character
        2. Raw sequence letters
        3. Quality values for the sequence

    In single cell data, the raw sequence often contains unique cell identifiers
    and unique molecular identifiers (UMI). These identifiers can the trimmed
    from the raw sequence.
    """

    def __init__(self, ident, seq, qual):
        """Initialization of a single read.

        Args:
            ident (str): read identifier.
            seq (str): read sequence.
            qual (str): read quality values.

        returns:
            None.
        """
        self.id = ident
        self.seq = seq.upper()
        self.qual = qual

        self.UMI = ''
        self.cell_ID = ''

    def __repr__(self):
        return self.id

    def __str__(self):
        return '{}\n{}\n+\n{}\n'.format(self.id, self.seq, self.qual)

    def __lt__(self, other):
        if self.id < other.id:
            return True
        return False

    def __eq__(self, other):
        if self.id == other.id:
            return True
        return False

    def set_UMI(self, UMI):
        """Sets the reads UMI
        Args:
            UMI (str): corresponding UMI.

        Returns:
            None.

        """
        self.UMI = UMI

    def set_cell_ID(self, cell_ID):
        """Sets the reads cell ID
        Args:
            cell_ID (str): corresponding cell ID.

        Returns:
            None.

        """
        self.cell_ID = cell_ID

    def trim_UMI(self, no_bp):
        """Trims the UMI from the sequence
        Args:
            no_bp (int): Number of base pairs of the UMI.

        Returns:
            str: UMI.

        """
        self.UMI = self.seq[:no_bp]
        self.seq = self.seq[no_bp:]
        return self.UMI

    def trim_cell_ID(self, no_bp):
        """Trims the cell ID from the sequence
        Args:
            no_bp (int): Number of base pairs of the cell ID.

        Returns:
            str: cell ID.

        """
        self.cell_ID = self.seq[:no_bp]
        self.seq = self.seq[no_bp:]
        return self.cell_ID

    def get_UMI(self):
        """Get the reads UMI

        Returns:
            str: UMI.

        """
        return self.UMI

    def get_cell_ID(self):
        """Get the reads cell ID

        Returns:
            str: cell ID.

        """
        return self.cell_ID

    def get_id(self):
        """Get the reads FASTQ id

        Returns:
            str: id.

        """
        return self.id

    def get_seq(self):
        """Get the reads Sequence

        Returns:
            str: Sequence.

        """
        return self.seq

    def get_qual(self):
        """Get the reads quality

        Returns:
            str: fastq quality.

        """
        return self.qual

# ------------------------------------------------------------------------------

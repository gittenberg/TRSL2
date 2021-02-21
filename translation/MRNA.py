import logging as log

# read model parameters
from translation.parameters import cr

"""
mRNA class definition

This module generates an mRNA object representing a generic (not sequence-specific) mRNA molecule.
It acts as the base class for the sequence-specific mRNA_spec.

Parameters:
      index - the ID of this particular mRNA object in the cell
      length - the length in nucleotides
      geneID - the ID of the gene belonging to this particular mRNA object (there might be more than one mRNA sharing the same geneID)

Test comment for review branch
"""


class MRNA:
    """
    class representing a single polysome
    """
    def __init__(self, index, length, geneID=None, ribosomes={}):
        """
        initializes one mRNA molecule
        """
        self.index = index          # counts the unique mRNA molecules; no biological meaning
        self.length = length        # length of mRNA in nts # http://bionumbers.hms.harvard.edu//bionumber.aspx?id=107678&ver=1
        self.geneID = geneID        # corresponds to sequence and proteinID; there might be more than one mRNA with this geneID
        self.ribosomes = ribosomes  # keys between 0 and self.length - 3, value None: no AA-tRNA, <value>: AA-tRNA of type <value>
        self.tic = False            # initiation time for a certain ribosome on this mRNA
        self.toc = False            # termination time for a certain ribosome on this mRNA
        self.tic_toc = []           # list of measured (tic, toc) pairs for this mRNA

    def detach_ribosome(self, pos):
        """
        detaches a ribosome without tRNA at position pos
        """
        if pos in self.ribosomes:
            del self.ribosomes[pos]
            # log.debug("detach_ribosome: successful at pos %s", pos)
            success = True
        else:
            self.ribosomes[pos] = None
            log.warning("detach_ribosome: unsuccessful at pos %s", pos)
            success = False
        return success

    def attach_ribosome_at_start(self):
        """
        attaches a ribosome without tRNA at position 0
        """
        if not self.ribosomes or min(self.ribosomes.keys()) > 3 * cr:
            self.ribosomes[0] = None
            # log.debug("attach_ribosome_at_start: successfully attached ribosome at pos 0")
            success = True
        else:
            log.debug("attach_ribosome_at_start: unsuccessful")
            success = False
        return success

    def first_position_occupied(self):
        """
        returns True iff the first 3 * cr nts of an mRNA are occupied by a ribosome
        """
        if not self.ribosomes:  # no ribosomes
            return False
        elif min(self.ribosomes.keys()) > 3 * cr:  # ribosomes behind position 30 nt
            # log.debug("first_position_occupied: ribosomes = %s", self.ribosomes)
            return False
        else:
            return True

    def next_range_free(self, pos, by=3):
        """
        returns True iff no ribosomes within the next by + 3 * cr positions from pos
        TODO: potentially time-consuming
        """
        return not any([ribo in range(pos + 1, pos + 1 + by + 3 * cr) for ribo in self.ribosomes])

    def find_max_free_range(self, pos):
        """
        returns maximum free range downstream from pos
        this range does not include the ribosome footprint
        """
        downstream_ribosomes = [ribo for ribo in self.ribosomes if ribo > pos]
        if downstream_ribosomes:
            next_ribo_pos = min(downstream_ribosomes)
            # log.debug("find_max_free_range: next_ribo_pos, pos = %s, %s", next_ribo_pos, pos)
            max_free_range = next_ribo_pos - pos
        else:
            max_free_range = 3 * cr + 3  # effectively infinite if no downstream ribosomes due to 3 ' utr
        return max_free_range

    def termination_condition(self):
        """
        returns True iff a ribosome is at the end of the mRNA.
        the ribosome footprint cr is ignored in this now.
        """
        if self.ribosomes:
            if max(self.ribosomes.keys()) + 3 >= self.length:
                return True
            else:
                return False
        else:
            return False

    def prev_range_free(self, pos, by=3):
        """
        returns True iff no ribosomes within the previous by+3*cr positions from pos
        TODO: potentially time-consuming
        """
        return not any([ribo in range(pos - 1, pos - 1 - by - 3 * cr, -1) for ribo in self.ribosomes])

    def translocate_ribosome(self, pos, by=3):
        """
        moves a ribosome at position pos by by positions 3'wards
        returns True if move was successful
        """
        if pos in self.ribosomes:
            # we assume the next 3*by positions are free to translocate
            self.ribosomes[pos + by] = self.ribosomes[pos]
            del self.ribosomes[pos]
            success = True
            # log.debug("translocate_ribosome: successfully moved ribosome to position %s", pos+by)
        else:
            log.warning("translocate_ribosome: cannot move ribosome: ribosome not found at %s", pos)
            success = False
        return success

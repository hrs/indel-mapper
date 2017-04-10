from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class Reference(object):

    def __init__(self, name, n20, sequence, pam, reads, padding=3):
        self.name = name
        self.n20 = n20.upper()
        self.sequence = sequence.upper()
        self.pam = pam.upper()
        self.reverse_complement_n20 = self._get_reverse_complement_n20()
        self.reads = reads
        self.reads_with_indels_near_the_cutsite = []
        self.padding = padding

        # save computation
        self.is_valid = self._compute_is_valid()
        if self.is_valid:
            self.reads_with_indels_near_the_cutsite = self._compute_reads_with_indels_near_the_cutsite(self.padding)

    def is_ngg(self):
        return self.pam == "NGG"

    def _get_reverse_complement_n20(self):
        return str(Seq(self.n20, generic_dna).reverse_complement())

    def _compute_is_valid(self):
        # is n20 or reverse complement of n20 in sequence
        if self.n20 in self.sequence:
            low_index, high_index = self._get_low_high_n20_index(self.n20)

        elif self.reverse_complement_n20 in self.sequence:
            low_index, high_index = self._get_low_high_n20_index(self.reverse_complement_n20)

        else:
            return False

        if self.is_ngg():
            return low_index != -1 and high_index != -1 and self.sequence[high_index+2:high_index+4] == "GG"
        else:
            return low_index != -1 and high_index != -1 and self.sequence[low_index-3:low_index-1] == "CC"

    def _get_low_high_n20_index(self, n20_sequence):
        low_index = self.sequence.find(n20_sequence)
        high_index = low_index + len(n20_sequence) - 1
        return low_index, high_index

    def n20_pam_index(self):
        # finds the index between the n20 and the pam
        if self.n20 in self.sequence:
            low_index, high_index = self._get_low_high_n20_index(self.n20)

        elif self.reverse_complement_n20 in self.sequence:
            low_index, high_index = self._get_low_high_n20_index(self.reverse_complement_n20)

        if self.is_ngg():
            return high_index
        else:
            return low_index

    def cutsite_index(self):
        if self.is_ngg():
            return self.n20_pam_index() - 3
        else:
            return self.n20_pam_index() + 3

    def pam_index(self):
        if self.is_ngg():
            return self.n20_pam_index() + 3
        else:
            return self.n20_pam_index() - 3

    def n20_index(self):
        if self.is_ngg():
            return self.n20_pam_index() - len(self.n20)
        else:
            return self.n20_pam_index() + len(self.n20)

    def _min_indel_dist(self, read):
        return min([self.distance_to_cutsite(i) for i in read.indels])

    def distance_to_cutsite(self, indel):

        if indel.is_enveloping_cutsite(self.cutsite_index()):
            return 0

        if indel.end_index < self.cutsite_index():
            return indel.end_index - self.cutsite_index()
        else:
            return indel.start_index - self.cutsite_index()

    def _min_abs_indel_dist(self, read):
        return min([abs(self.distance_to_cutsite(i)) for i in read.valid_indels])

    def _compute_reads_with_indels_near_the_cutsite(self, padding):
        # returns the reads with indels near the cutsite
        reads_with_indels = [read for read in self.reads if len(read.valid_indels) > 0]
        if len(reads_with_indels) > 0:
            return [read for read in reads_with_indels if self._min_abs_indel_dist(read) <= padding]
        else:
            return []

    def has_reads_with_indels_near_the_cutsite(self):
        return len(self.reads_with_indels_near_the_cutsite) > 0

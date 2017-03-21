from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class Reference(object):

    def __init__(self, name, n20, sequence, pam, reads):
        self.name = name
        self.n20 = n20.upper()
        self.sequence = sequence.upper()
        self.pam = pam.upper()
        self.reverse_complement_n20 = self._get_reverse_complement_n20()
        self.reads = reads

    def is_ngg(self):
        return self.pam == "NGG"

    def _get_reverse_complement_n20(self):
        return str(Seq(self.n20, generic_dna).reverse_complement())

    def is_n20_in_sequence(self):
        if self.n20 in self.sequence:
            low_index, high_index = self._get_low_high_n20_index(self.n20)

        if self.reverse_complement_n20 in self.sequence:
            low_index, high_index = self._get_low_high_n20_index(self.reverse_complement_n20)

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

        if self.reverse_complement_n20 in self.sequence:
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

    def distance_to_cutsite(self, indel):

        if indel.end_index < self.cutsite_index():
            if self.is_ngg():
                return self.cutsite_index() - indel.end_index
            else:
                return indel.end_index - self.cutsite_index()
        else:
            if self.is_ngg():
                return self.cutsite_index() - indel.start_index
            else:
                return indel.start_index - self.cutsite_index()

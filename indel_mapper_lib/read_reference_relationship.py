class ReadReferenceRelationship(object):
    def __init__(self, aligned_pair_index, aligned_pairs, reference_sequence, read_sequence, pam_index, n20_index, is_ngg):
        self.aligned_pair_index = aligned_pair_index
        self.aligned_pairs = aligned_pairs
        self.sequence_indexes = aligned_pairs[self.aligned_pair_index]
        self.read_index = self.sequence_indexes[0]
        self.reference_index = self.sequence_indexes[1]
        self.reference_sequence = reference_sequence
        self.read_sequence = read_sequence
        self.pam_index = pam_index
        self.n20_index = n20_index
        self.is_ngg = is_ngg

    def is_insertion(self):
        return self.reference_index is None

    def is_deletion(self):
        return self.read_index is None

    def _is_indel(self):
        return self.is_insertion() or self.is_deletion()

    def is_mismatch_but_not_indel(self):
        # this will error if it's an indel
        read_base = self.read_sequence[self.read_index]
        reference_base = self.reference_sequence[self.reference_index]
        return self._is_mismatch(read_base, reference_base)

    def is_between_pam_and_n20(self):
        if not self.is_ngg:
            return self.reference_index >= min([self.n20_index, self.pam_index]) and \
                self.reference_index < max([self.n20_index, self.pam_index])
        else:
            return self.reference_index >= min([self.n20_index+1, self.pam_index+1]) and \
                self.reference_index < max([self.n20_index+1, self.pam_index+1])

    def next_is_mismatch_or_indel(self):
        next_index = self.aligned_pair_index + 1
        if next_index < len(self.aligned_pairs):
            next_aligned_pair = self.aligned_pairs[next_index]
            return self._is_insertion(next_aligned_pair) or \
                self._is_deletion(next_aligned_pair) or \
                self._is_mismatch(self.read_sequence[next_aligned_pair[0]], self.reference_sequence[next_aligned_pair[1]])
        return False

    def _is_insertion(self, aligned_pair):
        return aligned_pair[1] is None

    def _is_deletion(self, aligned_pair):
        return aligned_pair[0] is None

    def _is_mismatch(self, read_base, reference_base):
        return read_base != reference_base

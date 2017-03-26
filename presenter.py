class SequenceTally(object):

    def __init__(self, read_sequence_presentation, reference_sequence_presentation, read_name):
        self.read_sequence_presentation = read_sequence_presentation
        self.reference_sequence_presentation = reference_sequence_presentation
        self.read_names = [read_name]

    def add_read(self, read_name):
        self.read_names.append(read_name)

    def _get_count(self):
        return len(self.read_names)

    def display_values(self):
        return ["count: {}".format(self._get_count()),
                "reads: {}".format(', '.join(self.read_names)),
                "ref : {}".format(self.reference_sequence_presentation),
                "read: {}".format(self.read_sequence_presentation)]

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
        return self._is_insertion(self.sequence_indexes)

    def is_deletion(self):
        return self._is_deletion(self.sequence_indexes)

    def _is_indel(self):
        return self.is_insertion() or self.is_deletion()

    def is_mismatch(self):
        if not self._is_indel():
            read_base = self.read_sequence[self.read_index]
            reference_base = self.reference_sequence[self.reference_index]
            return self._is_mismatch(read_base, reference_base)
        return True

    def is_between_pam_and_n20(self):
        if not self.is_ngg:
            return self.reference_index >= min([self.n20_index, self.pam_index]) and \
                self.reference_index < max([self.n20_index, self.pam_index])
        else:
            return self.reference_index >= min([self.n20_index+1, self.pam_index+1]) and \
                self.reference_index < max([self.n20_index+1, self.pam_index+1])

    def next_to_mismatch_or_indel(self):
        return self._previous_is_indel_or_mismatch() or self._next_is_indel_or_mismatch()

    def _previous_is_indel_or_mismatch(self):
        prev_index = self.aligned_pair_index - 1
        if prev_index >= 0:
            prev_aligned_pair = self.aligned_pairs[prev_index]
            return self._is_insertion(prev_aligned_pair) or \
                self._is_deletion(prev_aligned_pair) or \
                self._is_mismatch(self.read_sequence[prev_aligned_pair[0]], self.reference_sequence[prev_aligned_pair[1]])
        return False

    def _next_is_indel_or_mismatch(self):
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

class Presenter(object):

    def __init__(self, references):
        self.references = references

    def present(self):
        return '\n'.join([self.display_indels_near_cutsite(reference) for reference in self.references if len(reference.reads_with_indels_near_the_cutsite()) > 0])

    def _present_name(self, reference):
        return "name: {}".format(reference.name)

    def _present_reference(self, reference):
        return "sequence: {}".format(reference.sequence)

    def _present_n20(self, reference):
        return "n20: {}".format(reference.n20)

    def _present_pam(self, reference):
        return "pam: {}".format(reference.pam)

    def display_indels_near_cutsite(self, reference):
        return '\n'.join([self._present_name(reference),
                          self._present_reference(reference),
                          self._present_n20(reference),
                          self._present_pam(reference)] + self._tallies_as_display_values(reference))

    def _tally_sequences(self, reference):
        tallies = {}
        for read in reference.reads_with_indels_near_the_cutsite():
            reference_presentation, read_presentation = self.present_sequence(reference, read)
            if read_presentation in tallies:
                tallies[read_presentation].add_read(read.query_name)
            else:
                tallies[read_presentation] = SequenceTally(read_presentation, reference_presentation, read.query_name)
        return tallies

    def _tallies_as_display_values(self, reference):
        tallies = self._tally_sequences(reference)
        display_values = []
        for read_presentation, tally in tallies.items():
            display_values += tally.display_values()
        return display_values

    def present_sequence(self, reference, read):

        reference_presentation, read_presentation = self.get_sequence_representation(reference, read)
        reference_presentation_with_sites, read_presentation_with_sites = self.denote_cas9_sites(
            reference_presentation, read_presentation, reference, read)
        return reference_presentation_with_sites, read_presentation_with_sites

    def get_sequence_representation(self, reference, read):
        aligned_pairs = read.aligned_pairs
        reference_sequence = reference.sequence
        read_sequence = read.query_sequence

        reference_presentation = []
        read_presentation = []

        match_marker = "-"
        indel_marker = "_"

        for aligned_pair_index, sequence_indexes in enumerate(aligned_pairs):
            read_index, reference_index = sequence_indexes

            relationship = ReadReferenceRelationship(aligned_pair_index,
                                                     aligned_pairs,
                                                     reference_sequence,
                                                     read_sequence,
                                                     reference.pam_index(),
                                                     reference.n20_index(),
                                                     reference.is_ngg())

            if relationship.is_insertion():
                reference_presentation.append(indel_marker)
                read_presentation.append(read_sequence[read_index])
            elif relationship.is_deletion():
                reference_presentation.append(reference_sequence[reference_index])
                read_presentation.append(indel_marker)
            else:
                read_base = read_sequence[read_index]
                reference_base = reference_sequence[reference_index]
                if relationship.is_mismatch() or \
                   relationship.is_between_pam_and_n20() or \
                   relationship.next_to_mismatch_or_indel():
                    reference_presentation.append(reference_base)
                    read_presentation.append(read_base)
                else:
                    reference_presentation.append(match_marker)
                    read_presentation.append(match_marker)
        return reference_presentation, read_presentation

    def denote_cas9_sites(self, reference_presentation, read_presentation, reference, read):
        # denotes the positions of the cutsite, the n20, and the pam

        modified_reference_presentation = reference_presentation
        modified_read_presentation = read_presentation

        cutsite_index = reference.cutsite_index()
        pam_index = reference.pam_index()
        n20_pam_index = reference.n20_pam_index()
        n20_index = reference.n20_index()

        present_cutsite_index, present_pam_index, present_n20_pam_index, present_n20_index = self.get_presentation_indexes(
            read.aligned_pairs, cutsite_index, pam_index, n20_pam_index, n20_index)

        border_marker = "|"
        cutsite_marker = "||"

        borders =  [present_pam_index, present_n20_pam_index, present_n20_index]
        filtered_borders = [border for border in borders if border is not None]

        if not reference.is_ngg():
            for each_index in filtered_borders:
                modified_reference_presentation.insert(each_index, border_marker)
                modified_read_presentation.insert(each_index, border_marker)
            if present_cutsite_index is not None:
                modified_reference_presentation.insert(present_cutsite_index, cutsite_marker)
                modified_read_presentation.insert(present_cutsite_index, cutsite_marker)
        else:
            for each_index in filtered_borders:
                modified_reference_presentation.insert(each_index + 1, border_marker)
                modified_read_presentation.insert(each_index + 1, border_marker)
            if present_cutsite_index is not None:
                modified_reference_presentation.insert(present_cutsite_index + 1, cutsite_marker)
                modified_read_presentation.insert(present_cutsite_index + 1, cutsite_marker)

        return ''.join(modified_reference_presentation), ''.join(modified_read_presentation)

    def get_presentation_indexes(self, aligned_pairs, cutsite_index, pam_index, n20_pam_index, n20_index):
        presentation_cutsite_index = None
        presentation_pam_index = None
        presentation_n20_pam_index = None
        presentation_n20_index = None
        for aligned_pair_index, sequence_indexes in enumerate(aligned_pairs):
            _, reference_index = sequence_indexes
            if reference_index == cutsite_index:
                presentation_cutsite_index = aligned_pair_index
            elif reference_index == pam_index:
                presentation_pam_index = aligned_pair_index
            elif reference_index == n20_pam_index:
                presentation_n20_pam_index = aligned_pair_index
            elif reference_index == n20_index:
                presentation_n20_index = aligned_pair_index

        return presentation_cutsite_index, presentation_pam_index, presentation_n20_pam_index, presentation_n20_index

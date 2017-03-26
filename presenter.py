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

            if self._is_insertion(sequence_indexes):
                reference_presentation.append(indel_marker)
                read_presentation.append(read_sequence[read_index])
            elif self._is_deletion(sequence_indexes):
                reference_presentation.append(reference_sequence[reference_index])
                read_presentation.append(indel_marker)
            else:
                read_base = read_sequence[read_index]
                reference_base = reference_sequence[reference_index]
                if self._is_mismatch(read_base, reference_base) or \
                   self._is_between_pam_and_n20(reference, reference_index) or \
                   self._next_to_mismatch_or_indel(aligned_pair_index, aligned_pairs, reference_sequence, read_sequence):
                    reference_presentation.append(reference_base)
                    read_presentation.append(read_base)
                else:
                    reference_presentation.append(match_marker)
                    read_presentation.append(match_marker)
        return reference_presentation, read_presentation

    def _next_to_mismatch_or_indel(self, aligned_pair_index, aligned_pairs, reference_sequence, read_sequence):
        return self._previous_is_indel_or_mismatch(aligned_pair_index, aligned_pairs, reference_sequence, read_sequence) or \
            self._next_is_indel_or_mismatch(aligned_pair_index, aligned_pairs, reference_sequence, read_sequence)

    def _previous_is_indel_or_mismatch(self, aligned_pair_index, aligned_pairs, reference_sequence, read_sequence):
        prev_index = aligned_pair_index - 1
        if prev_index >= 0:
            prev_aligned_pairs = aligned_pairs[prev_index]
            return self._is_insertion(prev_aligned_pairs) or \
                self._is_deletion(prev_aligned_pairs) or \
                self._is_mismatch(read_sequence[prev_aligned_pairs[0]], reference_sequence[prev_aligned_pairs[1]])
        return False

    def _next_is_indel_or_mismatch(self, aligned_pair_index, aligned_pairs, reference_sequence, read_sequence):
        next_index = aligned_pair_index - 1
        if next_index < len(aligned_pairs):
            next_aligned_pairs = aligned_pairs[next_index]
            return self._is_insertion(next_aligned_pairs) or \
                self._is_deletion(next_aligned_pairs) or \
                self._is_mismatch(read_sequence[next_aligned_pairs[0]], reference_sequence[next_aligned_pairs[1]])
        return False

    def _is_between_pam_and_n20(self, reference, reference_index):
        n20_index = reference.n20_index()
        pam_index = reference.pam_index()
        if not reference.is_ngg():
            return reference_index >= min([n20_index, pam_index]) and reference_index < max([n20_index, pam_index])
        else:
            return reference_index >= min([n20_index+1, pam_index+1]) and reference_index < max([n20_index+1, pam_index+1])

    def _is_mismatch(self, read_base, reference_base):
        return read_base != reference_base

    def _is_insertion(self, sequence_indexes):
        return sequence_indexes[1] is None

    def _is_deletion(self, sequence_indexes):
        return sequence_indexes[0] is None

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

    def present_sequence2(self, reference, read):
        aligned_pairs = read.aligned_pairs
        reference_sequence = reference.sequence
        read_sequence = read.query_sequence

        reference_presentation = ''
        read_presentation = ''

        cutsite_index = reference.cutsite_index()
        pam_index = reference.pam_index()
        n20_pam_index = reference.n20_pam_index()
        n20_index = reference.n20_index()

        for aligned_pair_index, sequence_indexes in enumerate(aligned_pairs):
            read_index, reference_index = sequence_indexes

            if not reference.is_ngg():
                if reference_index == cutsite_index:
                    read_presentation += "||"
                    reference_presentation += "||"
                if reference_index == pam_index:
                    read_presentation += "|"
                    reference_presentation += "|"
                if reference_index == n20_pam_index:
                    read_presentation += "|"
                    reference_presentation += "|"
                if reference_index == n20_index:
                    read_presentation += "|"
                    reference_presentation += "|"

            if read_index is not None and reference_index is None:
                read_presentation += read_sequence[read_index]
                reference_presentation += '_'
            elif read_index is None and reference_index is not None:
                read_presentation += "_"
                reference_presentation += reference_sequence[reference_index]
            else:
                if self.previous_is_indel_and_current_is_not_indel(aligned_pairs, aligned_pair_index):
                    read_presentation += read_sequence[read_index]
                    reference_presentation += reference_sequence[reference_index]
                elif self.current_is_not_indel_and_next_is_indel(aligned_pairs, aligned_pair_index):
                    read_presentation += read_sequence[read_index]
                    reference_presentation += reference_sequence[reference_index]
                else:
                    if not reference.is_ngg() and reference_index >= min([n20_index, pam_index]) and reference_index < max([n20_index, pam_index]):
                        read_presentation += read_sequence[read_index]
                        reference_presentation += reference_sequence[reference_index]
                    elif reference.is_ngg() and reference_index >= min([n20_index+1, pam_index+1]) and reference_index < max([n20_index+1, pam_index+1]):
                        read_presentation += read_sequence[read_index]
                        reference_presentation += reference_sequence[reference_index]
                    else:
                        if read_sequence[read_index] != reference_sequence[reference_index]:
                            read_presentation += read_sequence[read_index]
                            reference_presentation += reference_sequence[reference_index]
                        else:
                            read_presentation += "-"
                            reference_presentation += "-"

            if reference.is_ngg():
                if reference_index == cutsite_index:
                    read_presentation += "||"
                    reference_presentation += "||"
                if reference_index == pam_index:
                    read_presentation += "|"
                    reference_presentation += "|"
                if reference_index == n20_pam_index:
                    read_presentation += "|"
                    reference_presentation += "|"
                if reference_index == n20_index:
                    read_presentation += "|"
                    reference_presentation += "|"

        return reference_presentation, read_presentation

    def previous_is_indel_and_current_is_not_indel(self, aligned_pairs, index):
        prev_index = index - 1
        if prev_index >= 0:
            return self.is_indel(aligned_pairs[prev_index]) and self.is_not_indel(aligned_pairs[index])
        return False

    def current_is_not_indel_and_next_is_indel(self, aligned_pairs, index):
        next_index = index + 1
        if next_index < len(aligned_pairs):
            return self.is_indel(aligned_pairs[next_index]) and self.is_not_indel(aligned_pairs[index])
        return False

    def is_indel(self, pair):
        read_index, reference_index = pair
        return (read_index is None) or (reference_index is None)

    def is_not_indel(self, pair):
        read_index, reference_index = pair
        return (read_index is not None) and (reference_index is not None)

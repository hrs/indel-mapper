import re


class SequenceTally(object):

    def __init__(self, cutsite_region_representation, sequence_relationship_presentation):
        self.cutsite_region_representation = cutsite_region_representation
        self.presentations = [sequence_relationship_presentation]

    def add_read(self, presentation):
        self.presentations.append(presentation)

    def count(self):
        return len(self.presentations)


class SequenceRelationshipPresentation(object):

    def __init__(self, read_sequence_presentation, reference_sequence_presentation, read_name):
        self.read_sequence_presentation = read_sequence_presentation
        self.reference_sequence_presentation = reference_sequence_presentation
        self.read_name = read_name


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


class DenotationIndex(object):

    def __init__(self, index, cutsite=False):
        self.index = index
        self.cutsite = cutsite
        self.representation = self._representation()

    def _representation(self):
        if self.cutsite:
            return "||"
        else:
            return "|"


class Cas9Denotations(object):

    def __init__(self, cutsite_index, pam_index, n20_pam_index, n20_index, aligned_pairs, is_ngg):
        self.cutsite_index = cutsite_index
        self.pam_index = pam_index
        self.n20_pam_index = n20_pam_index
        self.n20_index = n20_index
        self.aligned_pairs = aligned_pairs
        self.is_ngg = is_ngg
        self.denotation_indexes = self._get_denotation_indexes()
        self.filtered_denotation_indexes = [index_object for index_object in self.denotation_indexes if index_object.index is not None]

    def _get_denotation_indexes(self):
        denotation_indexes = []
        for aligned_pair_index, sequence_indexes in enumerate(self.aligned_pairs):
            _, reference_index = sequence_indexes
            if reference_index == self.cutsite_index:
                denotation_indexes.append(DenotationIndex(aligned_pair_index, True))
            elif reference_index == self.pam_index or reference_index == self.n20_pam_index or reference_index == self.n20_index:
                denotation_indexes.append(DenotationIndex(aligned_pair_index))

        return denotation_indexes

    def apply_to_presentation(self, reference_presentation_array, read_presentation_array):

        reference_presentation_string = ''
        read_presentation_string = ''
        if len(self.filtered_denotation_indexes) > 0:
            for index, value in enumerate(reference_presentation_array):
                if not self.is_ngg:
                    for denotation in self.filtered_denotation_indexes:
                        if index == denotation.index:
                            reference_presentation_string += denotation.representation
                            read_presentation_string += denotation.representation

                reference_presentation_string += value
                read_presentation_string += read_presentation_array[index]

                if self.is_ngg:
                    for denotation in self.filtered_denotation_indexes:
                        if index == denotation.index:
                            reference_presentation_string += denotation.representation
                            read_presentation_string += denotation.representation
        else:
            reference_presentation_string = ''.join(reference_presentation_array)
            read_presentation_string = ''.join(read_presentation_array)

        return reference_presentation_string, read_presentation_string


class ReferencePresenter(object):

    def __init__(self, reference):
        self.reference = reference
        self.tallies = sorted(self._tally_sequences(reference).values(), key=lambda tally: tally.count(), reverse=True)

    def name(self):
        return self.reference.name

    def sequence(self):
        return self.reference.sequence

    def n20(self):
        return self.reference.n20

    def pam(self):
        return self.reference.pam

    def _tally_sequences(self, reference):
        tallies = {}
        for read in reference.reads_with_indels_near_the_cutsite:
            reference_presentation, read_presentation, cutsite_region_presentation = self.present_sequence(reference, read)
            sequence_relationship_presentation = SequenceRelationshipPresentation(read_presentation, reference_presentation, read.query_name)
            if cutsite_region_presentation in tallies:
                tallies[cutsite_region_presentation].add_read(sequence_relationship_presentation)
            else:
                tallies[cutsite_region_presentation] = SequenceTally(cutsite_region_presentation, sequence_relationship_presentation)
        return tallies

    def present_sequence(self, reference, read):

        reference_presentation, read_presentation = self.get_sequence_representation(reference, read)
        reference_presentation_with_sites, read_presentation_with_sites = self.denote_cas9_sites(
            reference_presentation, read_presentation, reference, read)
        cas9_region_presentation = self.compute_cas9_presentation(read_presentation_with_sites)
        return reference_presentation_with_sites, read_presentation_with_sites, cas9_region_presentation

    def compute_cas9_presentation(self, read_presentation_string):
        areas_of_interest = re.split("[-]+", read_presentation_string)
        for area_of_interest in areas_of_interest:
            # there should be at least one | in one of the sections
            if "|" in area_of_interest:
                return area_of_interest
        return ""

    def get_sequence_representation(self, reference, read):
        aligned_pairs = read.aligned_pairs
        reference_sequence = reference.sequence
        read_sequence = read.query_sequence

        reference_presentation = []
        read_presentation = []

        match_marker = "-"
        indel_marker = "_"

        cutsite_region_presentation = ''

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

        cutsite_index = reference.cutsite_index()
        pam_index = reference.pam_index()
        n20_pam_index = reference.n20_pam_index()
        n20_index = reference.n20_index()

        denotations = Cas9Denotations(cutsite_index, pam_index, n20_pam_index, n20_index, read.aligned_pairs, reference.is_ngg())

        return denotations.apply_to_presentation(reference_presentation, read_presentation)


class Presenter(object):

    def __init__(self, references):
        self.references = references

    def present(self):
        return [ReferencePresenter(reference) for reference in self.references if reference.has_reads_with_indels_near_the_cutsite]

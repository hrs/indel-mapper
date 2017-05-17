import re
from .alignment import Alignment
from .cas9_indicator_inserter import Cas9IndicatorInserter
from .read_reference_relationship import ReadReferenceRelationship
from .realigner import Realigner
from .mutation_cluster import MutationCluster

CUTSITE_REPRESENTATION = "||" # string representation of where the cutsite is
PAM_AND_N20_REPRESENTATION = "|" # string representation of the boundaries of the pam and n20

class ReferenceProcessor(object):
    def compute_mutation_clusters(self, reference):
        all_clusters = self._cluster_reads_by_mutations_near_cutsite(reference).values()
        interesting_clusters = [cluster for cluster in all_clusters if cluster.has_mutations_adjacent_to_cutsite()]
        return sorted(interesting_clusters,
                      key=lambda cluster: cluster.count(), reverse=True)

    def _cluster_reads_by_mutations_near_cutsite(self, reference):
        """Reference -> {String: MutationCluster}"""
        clusters = {}
        for read in reference.reads_with_indels_near_the_cutsite:
            marked_sequence_alignment, marked_cas9_region = self._get_marked_alignments(reference, read)
            if (CUTSITE_REPRESENTATION in marked_cas9_region.read):
                if marked_cas9_region.read in clusters:
                    clusters[marked_cas9_region.read].add_read(marked_sequence_alignment)
                else:
                    clusters[marked_cas9_region.read] = MutationCluster(marked_sequence_alignment, marked_cas9_region)
        return clusters

    def _get_marked_alignments(self, reference, read):
        """Reference, Read -> Alignment, String"""
        reference_as_array, read_as_array = self.get_sequence_representation(reference, read)
        reference_string_with_cas9_markers, read_string_with_cas9_markers = self.denote_cas9_sites(
            reference_as_array, read_as_array, reference, read)
        alignment_with_markers = Alignment(reference_string_with_cas9_markers, read_string_with_cas9_markers)
        cas9_region = self._compute_cas9_region(alignment_with_markers)
        return self._realign(alignment_with_markers, cas9_region)

    def _compute_cas9_region(self, alignment):
        """Alignment -> Alignment"""
        cas9_region_in_read = self._get_cas9_region(alignment.read)
        cas9_region_in_reference = self._get_cas9_region(alignment.reference)
        return Alignment(cas9_region_in_reference, cas9_region_in_read)

    def _get_cas9_region(self, read):
        """String -> String"""
        areas_of_interest = re.split("[-]+", read)
        for area_of_interest in areas_of_interest:
            # there should be at least one | in one of the sections
            if "|" in area_of_interest:
                return area_of_interest
        return ""

    def _realign(self, alignment, cas9_region):
        """Alignment, String -> Alignment, Alignment"""
        new_cas9_region = self._get_new_cas9_regions(alignment, cas9_region)

        new_alignment = Alignment(alignment.reference.replace(cas9_region.reference, new_cas9_region.reference),
                                  alignment.read.replace(cas9_region.read, new_cas9_region.read))
        return new_alignment, new_cas9_region

    def _get_new_cas9_regions(self, alignment, cas9_region):
        """Alignment, Alignment -> Alignment"""
        return Realigner(cas9_region).align()

    def get_sequence_representation(self, reference, read):
        """Reference, Read -> [Char], [Char]"""
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
        """[Char], [Char], Reference, Read -> String, String"""
        # denotes the positions of the cutsite, the n20, and the pam

        cutsite_index = reference.cutsite_index()
        pam_index = reference.pam_index()
        n20_pam_index = reference.n20_pam_index()
        n20_index = reference.n20_index()

        inserter = Cas9IndicatorInserter(cutsite_index, pam_index, n20_pam_index, n20_index, read.aligned_pairs, reference.is_ngg())

        return inserter.insert_indicators(reference_presentation, read_presentation)


class ReferencePresenter(object):
    def __init__(self, reference, mutation_clusters):
        self.reference = reference
        self.mutation_clusters = mutation_clusters

    def name(self):
        return self.reference.name

    def sequence(self):
        return self.reference.sequence

    def n20(self):
        return self.reference.n20

    def pam(self):
        return self.reference.pam

    def has_mutation_clusters(self):
        return len(self.mutation_clusters) != 0

    def total_reads(self):
        return len(self.reference.reads)

    def csv_row_prefix_cells(self):
        return [self.name(), self.total_reads()]


class Presenter(object):
    def __init__(self, references):
        self.references = references

    def present(self):
        processor = ReferenceProcessor()
        results = []
        for reference in self.references:
            if reference.has_reads_with_indels_near_the_cutsite:
                mutation_clusters = processor.compute_mutation_clusters(reference)
                results.append(ReferencePresenter(reference, mutation_clusters))
        return results

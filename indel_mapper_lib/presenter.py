from .reference_processor import ReferenceProcessor

CUTSITE_REPRESENTATION = "||" # string representation of where the cutsite is
PAM_AND_N20_REPRESENTATION = "|" # string representation of the boundaries of the pam and n20

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

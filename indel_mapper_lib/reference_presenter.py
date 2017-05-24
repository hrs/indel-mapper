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

    def to_dict(self):
        return {
            "name": self.name(),
            "sequence": self.sequence(),
            "n20": self.n20(),
            "pam": self.pam(),
            "has_mutation_clusters": self.has_mutation_clusters(),
            "total_reads": self.total_reads(),
            "mutation_clusters": [cluster.to_dict() for cluster in self.mutation_clusters]
        }

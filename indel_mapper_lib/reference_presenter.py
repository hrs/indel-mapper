from .mutation_cluster import MutationCluster

class ReferencePresenter(object):
    def __init__(self, reference=None, mutation_clusters=None, from_dict=None):
        if from_dict is not None:
            self._init_from_dict(from_dict)
        else:
            self.mutation_clusters = mutation_clusters
            self._name = reference.name
            self._sequence = reference.sequence
            self._n20 = reference.n20
            self._pam = reference.pam
            self._total_reads = len(reference.reads)

    def name(self):
        return self._name

    def sequence(self):
        return self._sequence

    def n20(self):
        return self._n20

    def pam(self):
        return self._pam

    def has_mutation_clusters(self):
        return len(self.mutation_clusters) != 0

    def total_reads(self):
        return self._total_reads

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

    def _init_from_dict(self, reference_presenter_dict):
        self._name = reference_presenter_dict["name"]
        self._sequence = reference_presenter_dict["sequence"]
        self._pam = reference_presenter_dict["pam"]
        self._n20 = reference_presenter_dict["n20"]
        self._total_reads = reference_presenter_dict["total_reads"]
        self.mutation_clusters = [MutationCluster(from_dict=cluster_dict) for cluster_dict in reference_presenter_dict["mutation_clusters"]]

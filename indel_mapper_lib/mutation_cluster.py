from .sequence_diff import SequenceDiff
from .alignment import Alignment

class MutationCluster(object):
    def __init__(self, marked_alignment=None, marked_cas9_region=None, from_dict=None):
        print(from_dict)
        if from_dict is not None:
            self._init_from_dict(from_dict)
        else:
            self.cas9_region = marked_cas9_region
            self.alignments = [marked_alignment]
            self.description = SequenceDiff(marked_cas9_region.reference, marked_cas9_region.read).description()
            self._count = len(self.alignments)

    def add_read(self, representation):
        self.alignments.append(representation)
        self._count = len(self.alignments)

    def count(self):
        return self._count

    def csv_row(self):
        return [self.cas9_region.reference, self.cas9_region.read, self.description, self.count()]

    def has_mutations_adjacent_to_cutsite(self):
        if "||" in self.cas9_region.read:
            index_of_cutsite = self.cas9_region.read.index("||")
            return self._has_unmatched_bases_left_of_cutsite(index_of_cutsite) or \
                self._has_unmatched_bases_right_of_cutsite(index_of_cutsite)
        elif self.alignments[0].starts_with(self.cas9_region):
            return self.cas9_region.has_unmatched_bases_at(0)
        else:
            return self.cas9_region.has_unmatched_bases_at(-1)

    def _has_unmatched_bases_left_of_cutsite(self, cutsite):
        index = cutsite - 1
        return index >= 0 and self.cas9_region.has_unmatched_bases_at(index)

    def _has_unmatched_bases_right_of_cutsite(self, cutsite):
        index = cutsite + 1
        return index < len(self.cas9_region.read) and self.cas9_region.has_unmatched_bases_at(index)

    def to_dict(self):
        return {
            "cas9_region": self.cas9_region.to_dict(),
            "count": self.count(),
            "description": self.description,
            "alignments": [alignment.to_dict() for alignment in self.alignments]
        }

    def _init_from_dict(self, mutation_cluster_dict):
        self._count = mutation_cluster_dict["count"]
        self.description = mutation_cluster_dict["description"]
        self.cas9_region = Alignment(from_dict=mutation_cluster_dict["cas9_region"])
        self.alignments = [Alignment(from_dict=alignment_dict) for alignment_dict in mutation_cluster_dict["alignments"]]

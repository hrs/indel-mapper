import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Presenter(object):

    def __init__(self, references):
        self.references = references

    def present_reads_with_indels_near_cutsite(self):
        for reference in self.references:
            if len(reference.sorted_reads()) > 0:
                self.plot_sorted_reads(reference)

            #for read in sorted_reads:
            #    print(self.formatted_read(read, reference.cutsite))

    def formatted_read(self, read, cutsite):
        aligned_pairs = read.aligned_pairs
        read_str = " " * aligned_pairs[0][1]
        for query, reference in aligned_pairs:
            if reference == cutsite:
                read_str += "|"
            elif reference is None:
                read_str += "+"
            elif query is None:
                read_str += "-"
            else:
                read_str += "_"
        return read_str

    def get_values(self, reference):
        x_insertions = []
        y_insertions = []
        x_deletions = []
        y_deletions = []

        for idx, read in enumerate(reference.sorted_reads()[:50]):
            for indel in read.indels:
                if indel.is_insertion:
                    x_insertions.append(idx)
                    y_insertions.append(reference.distance_to_cutsite(indel))
                else:
                    x_deletions.append(idx)
                    y_deletions.append(reference.distance_to_cutsite(indel))

        return x_insertions, y_insertions, x_deletions, y_deletions

    def plot_sorted_reads(self, reference):
         x_insertions, y_insertions, x_deletions, y_deletions = self.get_values(reference)

         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.scatter(x_insertions, y_insertions, c='r')
         ax.scatter(x_deletions, y_deletions, c='b')

         ax.set_xlabel('Reads (Sorted by Min. Indel Distance from Cutsite)')
         ax.set_ylabel('Indel Distance from Cutsite')

         plt.savefig("results/" + reference.name + ".png")

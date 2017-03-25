#import matplotlib.pyplot as plt

class Presenter(object):

    def __init__(self, references):
        self.references = references

    def present_reads_with_indels(self):
        for reference in self.references:
            if len(reference.sorted_reads()) > 0:
                self.plot_sorted_reads(reference)

            #for read in sorted_reads:
            #    print(self.formatted_read(read, reference.cutsite))

    def present(self):
        for reference in self.references:
            reads = reference.reads_with_indels_near_the_cutsite()
            if len(reads) > 0 :
                self.print_formatted(reference, reads)

    def print_formatted(self, reference, reads):
        print(reference.name)
        print(reference.n20)
        print(reference.sequence)
        print(reference.pam)
        print("-------")
        for r in reads:
            print(r.query_name)
            reference_presentation, read_presentation = self.present_sequence(reference, r)
            print(reference_presentation)
            print(read_presentation)
        print("=========")

    def present_sequence(self, reference, read):
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

    def get_values(self, reference):
        x_insertions = []
        y_insertions = []
        x_deletions = []
        y_deletions = []

        for idx, read in enumerate(reference.sorted_reads()):
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

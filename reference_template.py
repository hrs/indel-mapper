class TemplateReference(object):

    def __init__(self, name, n20, sequence, pam):
        self.name = name
        self.n20 = n20
        self.sequence = sequence
        self.pam = pam

    def is_ngg(self):
        return self.pam == "NGG"

    def n20_index(self):
        # finds the index between the n20 and the pam
        if self.is_ngg(): # find the index of where n20 ends when NGG
            index = self.sequence.rfind(self.n20)
            if index == -1:
                n20_seq = Seq(self.n20, generic_dna)
                index = self.sequence.rfind(str(n20_seq.reverse_complement()))

        else:
            index = self.sequence.find(self.n20) # find index of where n20 begins when CCN
            if index == -1:
                n20_seq = Seq(self.n20, generic_dna)
                index = self.sequence.find(str(n20_seq.reverse_complement()))

        return index

    def cutsite_index(self):
        if self.is_ngg():
            return self.n20_index() - 3
        else:
            return self.n20_index() + 3

    def distance_to_cutsite(self, indel):

        if indel.end_index < self.cutsite_index():
            if self.is_ngg():
                return self.cutsite_index() - indel.end_index
            else:
                return indel.end_index - self.cutsite_index()
        else:
            if self.is_ngg():
                return self.cutsite_index() - indel.start_index
            else:
                return indel.start_index - self.cutsite_index()

from indel_mapper_lib.reference import Reference

class ReferenceParser():

    def __init__(self, reference_file, reads):
        self.reference_file = reference_file
        self.reads = reads

    def references(self):
        references_from_file = []

        for line in self.reference_file:
            name, n20, sequence, pam = line
            reads_for_reference = self._reads_for_reference(name)
            references_from_file.append(Reference(name.strip(), n20.strip(), sequence.strip(), pam.strip(), reads_for_reference))

        return references_from_file

    def _reads_for_reference(self, reference_name):
        return [read for read in self.reads if read.reference_name == reference_name]

from reference import Reference

class ReferenceParser():

    def __init__(self, reference_file):
        self.reference_file = reference_file

    def references(self):
        references_from_file = []

        for line in self.reference_file:
            windows_formatted_lines = line.split('\r')
            for wl in windows_formatted_lines:
                name, n20, sequence, pam = wl.split('\t')
                references_from_file.append(Reference(name, n20, sequence, pam))

        return references_from_file

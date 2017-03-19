# reads = SamParser(pysam.AlignmentFile(my_file_name, "rb")).reads()

class SamParser(object):

    def __init__(self, samfile):
        self.samfile = samfile

    def reads(self):
        return list(samfile.fetch())

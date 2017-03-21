# reads = SamParser(pysam.AlignmentFile(my_file_name, "rb")).reads()

from read import Read

class SamParser(object):

    def __init__(self, samfile):
        self.samfile = samfile

    def reads(self):
        return [self._make_read(filtered_read) for filtered_read in filter(self.filter_bad_reads, self._fetch())]

    def filter_bad_reads(self, pysam_read):
        try:
            pysam_read.reference_name # try and see if the reference name exists
            return len(pysam_read.get_reference_positions(full_length=True)) != 0
        except:
            return False

    def _fetch(self):
        return self.samfile.fetch()

    def _make_read(self, pysam_read):
        return Read(query_name=pysam_read.query_name,
                    reference_name=pysam_read.reference_name,
                    reference_positions=pysam_read.get_reference_positions(full_length=True))

    def _read_has_non_insertions_or_deletions(self, pysam_read):
        for cig in pysam_read.cigartuples:
            cig_type, _ = cig
            if self.is_other(cig_type):
                return True
        return False

    def is_insertion_or_deletion(self, cigar_type):
        return cigar_type == 1 or cigar_type == 2

    def is_other(self, cigar_type):
        return cigar_type != 0 and not self.is_insertion_or_deletion(cigar_type)

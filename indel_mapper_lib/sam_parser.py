from indel_mapper_lib.read import Read

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
                    reference_positions=pysam_read.get_reference_positions(full_length=True),
                    query_sequence=pysam_read.query_sequence,
                    aligned_pairs=pysam_read.get_aligned_pairs())

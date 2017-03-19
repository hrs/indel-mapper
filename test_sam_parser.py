import unittest
from sam_parser import SamParser

class FakeAlignmentFile(object):

    def __init__(self, reads):
        self.reads = reads

    def fetch(self):
        return self.reads

class TestSamParser(unittest.TestCase):

    def test_reads(self):
        expected_reads = [1,2,3]
        sam_file = FakeAlignmentFile(expected_reads)
        reads = SamParser(sam_file).reads()

        self.assertEqual(reads, expected_reads)

import unittest
from unittest.mock import patch
from sam_parser import SamParser

class TestSamParser(unittest.TestCase):

    def mock_fetch(self):
        class FakeAlignedSegment(object):
            def __init__(self, query_name=None, reference_name=None, reference_positions=None, cigartuples=()):
                self.query_name = query_name
                self.cigartuples = cigartuples
                if reference_name:
                    self.reference_name = reference_name
                if reference_positions:
                    self.reference_positions = reference_positions
            def get_reference_positions(self, full_length):
                return self.reference_positions

        segment_a = FakeAlignedSegment(query_name='foo')
        segment_b = FakeAlignedSegment(query_name='bar',
                                       reference_name='cats',
                                       reference_positions=())
        segment_c = FakeAlignedSegment(query_name='baz',
                                       reference_name='dogs',
                                       reference_positions=(1,2,None,None,3),
                                       cigartuples=((0,1),(0,2),(2,None),(2,None),(0,3)))

        return [segment_a, segment_b, segment_c]

    @patch('sam_parser.SamParser._fetch')
    def test_reads(self, mock_fetcher):
        sam_file = "fake alignment file"
        mock_values = self.mock_fetch()
        mock_fetcher.return_value = mock_values
        reads = SamParser(sam_file).reads()

        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].query_name, mock_values[-1].query_name)
        self.assertEqual(reads[0].reference_name, mock_values[-1].reference_name)
        self.assertEqual(reads[0].reference_positions, mock_values[-1].reference_positions)

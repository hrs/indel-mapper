import unittest
from reference_parser import ReferenceParser

class TestReferenceParser(unittest.TestCase):

    def test_reads(self):
        fake_windows_file = ['foo\tcat\tccaattgg\tngg\rbar\ttag\tggaattttg\tngg']
        references = ReferenceParser(fake_windows_file).references()

        self.assertEqual(len(references), 2)

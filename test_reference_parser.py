import unittest
from reference_parser import ReferenceParser
from read import Read

class TestReferenceParser(unittest.TestCase):

    def test_reads(self):

        reads = [Read(query_name='lion',reference_name='foo', reference_positions=[1,2]),
                 Read(query_name='wolf',reference_name='foo', reference_positions=[4,5])]

        fake_windows_file = ['foo\tcat\tccaattgg\tngg\rbar\ttag\tggaattttg\tngg']
        references = ReferenceParser(fake_windows_file, reads).references()

        self.assertEqual(len(references), 2)

        foo_reference = references[0]
        bar_reference = references[1]

        self.assertEqual(len(foo_reference.reads), 2)
        self.assertEqual(foo_reference.reads[0].query_name, reads[0].query_name)
        self.assertEqual(foo_reference.reads[1].query_name, reads[1].query_name)
        self.assertEqual(len(bar_reference.reads), 0)

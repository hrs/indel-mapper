import unittest
from read import Read

class TestRead(unittest.TestCase):

    def test_create_read(self):

        query_name = "foo"
        reference_name = "lion"
        reference_positions = [1,2,3,4,None,None,5]

        test_read = Read(query_name, reference_name, reference_positions)

        self.assertEqual(test_read.query_name, query_name)
        self.assertEqual(test_read.reference_name, reference_name)
        self.assertEqual(test_read.reference_positions, reference_positions)


    def test_create_indels_for_read(self):
        query_name = ""
        reference_name = ""
        reference_positions_a = [1,2,3,4,None,None,None,5,6,7,None,None,None,None,8]
        read_a = Read(query_name, reference_name, reference_positions_a)
        indels_a = read_a.indels

        reference_positions_b = [3,4,5,10,11,12,None,None,13]
        read_b = Read(query_name, reference_name, reference_positions_b)
        indels_b = read_b.indels

        self.assertEqual(len(indels_a), 2)

        self.assertEqual(indels_a[0].start_index, 4)
        self.assertEqual(indels_a[0].end_index, 5)
        self.assertEqual(indels_a[0].length, 3)
        self.assertEqual(indels_a[0].is_deletion, False),
        self.assertEqual(indels_a[0].is_insertion, True),

        self.assertEqual(indels_a[1].start_index, 7)
        self.assertEqual(indels_a[1].end_index, 8)
        self.assertEqual(indels_a[1].length, 4)
        self.assertEqual(indels_a[1].is_deletion, False),
        self.assertEqual(indels_a[1].is_insertion, True),

        self.assertEqual(len(indels_b), 2)

        self.assertEqual(indels_b[0].start_index, 5)
        self.assertEqual(indels_b[0].end_index, 10)
        self.assertEqual(indels_b[0].length, 4)
        self.assertEqual(indels_b[0].is_deletion, True),
        self.assertEqual(indels_b[0].is_insertion, False),

        self.assertEqual(indels_b[1].start_index, 12)
        self.assertEqual(indels_b[1].end_index, 13)
        self.assertEqual(indels_b[1].length, 2)
        self.assertEqual(indels_b[1].is_deletion, False),
        self.assertEqual(indels_b[1].is_insertion, True),

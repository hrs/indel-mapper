import unittest
from indel import Indel

class TestIndel(unittest.TestCase):

    def test_create_indel(self):
        start_index = 3
        end_index = 10
        length = 6
        is_deletion = False

        indel = Indel(start_index, end_index, length, is_deletion)

        self.assertEqual(indel.start_index, start_index)
        self.assertEqual(indel.end_index, end_index)
        self.assertEqual(indel.length, length)
        self.assertFalse(indel.is_deletion)
        self.assertTrue(indel.is_insertion)

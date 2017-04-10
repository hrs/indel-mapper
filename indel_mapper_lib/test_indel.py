import unittest
from .indel import Indel

class TestIndel(unittest.TestCase):

    def test_create_indel(self):
        start_index = 3
        end_index = 10
        length = 6
        is_deletion = True

        indel = Indel(start_index, end_index, length, is_deletion)

        self.assertEqual(indel.start_index, start_index)
        self.assertEqual(indel.end_index, end_index)
        self.assertEqual(indel.length, length)
        self.assertTrue(indel.is_deletion)
        self.assertFalse(indel.is_insertion)

    def test_is_enveloping_cutsite(self):
        start_index = 3
        end_index = 10
        length = 6
        is_deletion = True

        indel = Indel(start_index, end_index, length, is_deletion)

        self.assertTrue(indel.is_enveloping_cutsite(4))
        self.assertFalse(indel.is_enveloping_cutsite(2))
        self.assertFalse(indel.is_enveloping_cutsite(11))

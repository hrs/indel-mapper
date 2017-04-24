import unittest

from .diff_components import DiffDeletion
from .diff_components import DiffInsertion
from .diff_components import DiffMatch
from .diff_components import DiffMetadata


class TestDiffDeletion(unittest.TestCase):
    def test_is_match(self):
        self.assertFalse(DiffDeletion("ACT").is_match())

    def test_is_metadata(self):
        self.assertFalse(DiffDeletion("ACT").is_metadata())

    def test_repr(self):
        self.assertEqual(str(DiffDeletion("A")),
                         "1 deletion (A)")
        self.assertEqual(str(DiffDeletion("ACT")),
                         "3 deletions (ACT)")


class TestDiffInsertion(unittest.TestCase):
    def test_is_match(self):
        self.assertFalse(DiffInsertion("ACT").is_match())

    def test_is_metadata(self):
        self.assertFalse(DiffInsertion("ACT").is_metadata())

    def test_repr(self):
        self.assertEqual(str(DiffInsertion("A")),
                         "1 insertion (A)")
        self.assertEqual(str(DiffInsertion("ACT")),
                         "3 insertions (ACT)")


class TestDiffMatch(unittest.TestCase):
    def test_is_match(self):
        self.assertTrue(DiffMatch("ACT").is_match())

    def test_is_metadata(self):
        self.assertFalse(DiffMatch("ACT").is_metadata())

    def test_repr(self):
        self.assertEqual(str(DiffMatch("A")), "1 match")
        self.assertEqual(str(DiffMatch("ACT")), "3 matches")


class TestDiffMetadata(unittest.TestCase):
    def test_is_match(self):
        self.assertFalse(DiffMetadata("||").is_match())

    def test_is_metadata(self):
        self.assertTrue(DiffMetadata("||").is_metadata())

    def test_repr(self):
        self.assertEqual(str(DiffMetadata("||")), "||")

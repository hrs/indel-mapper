import unittest

from .sequence_diff import SequenceDiff
from .diff_components import DiffDeletion
from .diff_components import DiffInsertion
from .diff_components import DiffMatch
from .diff_components import DiffMetadata


class TestSequenceDiff(unittest.TestCase):
    def test_validation(self):
        with self.assertRaises(ValueError):
            SequenceDiff("AT", "TAT")
        with self.assertRaises(ValueError):
            SequenceDiff("A||TA", "AT||A")

    def test_tokenizer(self):
        diff = SequenceDiff("GATT___ACA||AGGT__C",
                            "GAT_GGCA__||AGG_AGC")
        expected_changes = [DiffMatch("GAT"),
                            DiffDeletion("T"),
                            DiffInsertion("GGC"),
                            DiffMatch("A"),
                            DiffDeletion("CA"),
                            DiffMetadata("||"),
                            DiffMatch("AGG"),
                            DiffDeletion("T"),
                            DiffInsertion("AG"),
                            DiffMatch("C")]
        self.assertEqual(diff.changes, expected_changes)

    def test_description(self):
        diff = SequenceDiff("GATT___ACA||AGGT__C",
                            "GAT_GGCA__||AGG_AGC")
        self.assertEqual(diff.description(),
                         "1 deletion (T), 3 insertions (GGC), 1 match, "
                         "2 deletions (CA) || 3 matches, 1 deletion (T), "
                         "2 insertions (AG)")

    def test_description_without_cutsite(self):
        diff = SequenceDiff("GATT___ACAAGGT__C",
                            "GAT_GGCA__AGG_AGC")
        self.assertEqual(diff.description(),
                         "1 deletion (T), 3 insertions (GGC), 1 match, "
                         "2 deletions (CA), 3 matches, 1 deletion (T), "
                         "2 insertions (AG)")

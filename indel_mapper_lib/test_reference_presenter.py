import unittest
from .reference import Reference
from .reference_presenter import ReferencePresenter

class TestReferencePresenter(unittest.TestCase):

    def test_reference_presenter(self):
        ref = ReferencePresenter(Reference("foo", "", "", "", []), [])

        self.assertEqual(ref.csv_row_prefix_cells(), ["foo", 0])

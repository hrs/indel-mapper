import unittest
from .alignment import Alignment

class TestAlignment(unittest.TestCase):
    def test_equality(self):
        self.assertEqual(Alignment("foo", "bar"),
                         Alignment("foo", "bar"))
        self.assertNotEqual(Alignment("foo", "bar"),
                            Alignment("bar", "foo"))

    def test_repr(self):
        self.assertEqual(str(Alignment("foo", "bar")),
                         "(foo, bar)")

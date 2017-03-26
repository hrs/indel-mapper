import unittest
from read import Read
from presenter import SequenceTally
from presenter import Presenter

class TestSequenceTally(unittest.TestCase):

    def test_sequence_tally(self):

        reference_presentation = "aaaa"
        read_presentation = "aaat"
        read_name = "foo"

        tally = SequenceTally(read_presentation, reference_presentation, read_name)

        tally.add_read("bar")
        tally.add_read("baz")

        display_values = tally.display_values()

        self.assertEqual(len(display_values), 4)
        self.assertEqual(display_values[0], "count: 3")
        self.assertEqual(display_values[1], "reads: foo, bar, baz")
        self.assertEqual(display_values[2], "ref : {}".format(reference_presentation))
        self.assertEqual(display_values[3], "read: {}".format(read_presentation))

class TestPresenter(unittest.TestCase):

    def test_foo(self):
        pass

import unittest
from reference import Reference

class TestReference(unittest.TestCase):

    def test_create_reference(self):
        name = "foo"
        n20 = "aaaaggg"
        sequence = "catcatcatcat"
        pam = "ngg"

        reference = Reference(name, n20, sequence, pam)

        self.assertEqual(reference.name, name)
        self.assertEqual(reference.n20, n20.upper())
        self.assertEqual(reference.sequence, sequence.upper())
        self.assertEqual(reference.pam, pam.upper())
        self.assertEqual(reference.reverse_complement_n20, "CCCTTTT")

    def test_is_ngg_pam(self):
        ngg_reference = Reference("", "", "", "ngg")
        ccn_reference = Reference("", "", "", "ccn")

        self.assertTrue(ngg_reference.is_ngg())
        self.assertFalse(ccn_reference.is_ngg())

    def test_is_n20_in_sequence(self):
        ngg_reference = Reference("", "cat", "ccattgg", "ngg")
        bad_ngg_reference = Reference("", "cat", "ccatt", "ngg")
        ccn_reference = Reference("", "cat", "ccgcatt", "ccg")
        bad_ccn_reference = Reference("", "cat", "ttccatt", "ccg")
        rc_ngg_reference = Reference("", "atg", "ccattgg", "ngg")
        rc_ccn_reference = Reference("", "atg", "ccgcatt", "cgg")

        self.assertTrue(ngg_reference.is_n20_in_sequence())
        self.assertFalse(bad_ngg_reference.is_n20_in_sequence())
        self.assertTrue(ccn_reference.is_n20_in_sequence())
        self.assertFalse(bad_ccn_reference.is_n20_in_sequence())
        self.assertTrue(rc_ngg_reference.is_n20_in_sequence())
        self.assertTrue(rc_ccn_reference.is_n20_in_sequence())

    def test_cutsite(self):
        n20 = "aaaatttc"
        sequence = "tactactacaaaatttcnggt"
        pam = "ngg"

        ngg_reference = Reference("", n20, sequence, pam)

        n20 = "cattttg"
        sequence = "ccgcattttgaaa"
        pam = "ccn"

        ccn_reference = Reference("", n20, sequence, pam)

        self.assertEqual(ngg_reference.n20_pam_index(), 16)
        self.assertEqual(ngg_reference.cutsite_index(), 13)
        self.assertEqual(ccn_reference.n20_pam_index(), 3)
        self.assertEqual(ccn_reference.cutsite_index(), 6)

import unittest
from read import Read
from presenter import SequenceTally
from presenter import Presenter
from presenter import ReadReferenceRelationship
from presenter import DenotationIndex
from presenter import Cas9Denotations

class TestSequenceTally(unittest.TestCase):

    def test_sequence_tally(self):

        reference_presentation = "aaaa"
        read_presentation = "aaat"
        read_name = "foo"

        tally = SequenceTally(read_presentation, reference_presentation, read_name)

        tally.add_read("bar")
        tally.add_read("baz")

        self.assertEqual(tally.read_names_as_string(), "foo, bar, baz")
        self.assertEqual(tally.count(), 3)
        self.assertEqual(tally.reference_sequence_presentation, reference_presentation)
        self.assertEqual(tally.read_sequence_presentation, read_presentation)

class TestReadReferenceRelationship(unittest.TestCase):

    def create_relationship(self, aligned_pair_index, is_ngg=True):

        aligned_pairs = [(0, 14), (1, 15), (2, 16), (3, None), (4, None), (5, None), (6, 17), (7, 18), (8, 19), (9, 20), (10, 21), (11, 22), (12, 23), (13, 24), (14, 25), (15, 26), (16, 27), (17, 28), (18, 29), (19, 30), (20, 31), (21, 32), (22, 33), (23, 34), (24, 35), (25, 36), (26, 37), (27, 38), (28, 39), (29, 40), (30, 41), (31, 42), (32, 43), (33, 44), (34, 45), (35, 46), (36, 47), (37, 48), (38, 49), (39, 50), (40, 51), (41, 52), (42, 53), (43, 54), (44, 55), (45, 56), (46, 57), (47, 58), (48, 59), (49, 60), (50, 61), (51, 62), (52, 63), (53, 64), (54, 65), (55, 66), (56, 67), (57, 68), (58, 69), (59, 70), (60, 71), (61, 72), (62, 73), (63, 74), (64, 75), (65, 76), (66, 77), (67, 78), (68, 79), (69, 80), (70, 81), (71, 82), (72, 83), (73, 84), (74, 85), (75, 86), (76, 87), (77, 88), (78, 89), (79, 90), (80, 91), (81, 92), (82, 93), (83, 94), (84, 95), (85, 96), (86, 97), (87, 98), (88, 99), (89, 100), (90, 101), (91, 102), (92, 103), (93, 104), (94, 105), (95, 106), (96, 107), (97, 108), (98, 109), (99, 110), (100, 111), (101, 112), (102, 113), (103, 114), (104, 115), (105, 116), (None, 117), (106, 118), (107, 119), (108, 120), (109, 121), (110, 122), (111, 123), (112, 124), (113, 125), (114, 126), (115, 127), (116, 128), (117, 129), (118, 130), (119, 131)]
        reference_sequence = "AGGAACTGGGAGAGGACGATCCGGTTAGGGAGGTTGGGGAACTAATCTCAACGCTGCGTTTACAGATGAAGCCGCTTTTATATGGCGTATATGTTTGCTTAGAGGGGCCGACGGAGATTAGGAGAAGCCATCCTTTGGCGCCAATGATCAAAGCGTCTGCCAAGGAGAAGAAGCCAAGGGATGGGCCTTTCAGAGAGGGCAAGGAGTCATGCTGCTCTGGATGCCAGTGTCAGGACAAG"
        read_sequence = "TCCGCAGATCCGGTTAGGGAGGTTGGGGAACTAATCTCAACGCTGCGTTTACAGATGAAGCCGCTTTTATATGGCGTATATGTTTGCTTAGAGGGGCCGACGGAGATAGGAGAAGCCATCCTTTGGCGCCAATGATCAAAGCGTCTGCCAAGGAGAAGAAGCCAAGGGATGGGCCTTTCAGAGAGGGCAAGGAGTCATGCTGCTCTGGATGCCAGTGTCAGGACAA"
        pam_index = 121
        n20_index = 98

        return ReadReferenceRelationship(aligned_pair_index,
                                         aligned_pairs,
                                         reference_sequence,
                                         read_sequence,
                                         pam_index,
                                         n20_index,
                                         is_ngg)

    def test_is_insertion(self):
        relationship = self.create_relationship(aligned_pair_index=4)
        self.assertTrue(relationship.is_insertion())

        relationship = self.create_relationship(aligned_pair_index=2)
        self.assertFalse(relationship.is_insertion())

        relationship = self.create_relationship(aligned_pair_index=106)
        self.assertFalse(relationship.is_insertion())

    def test_is_deletion(self):
        relationship = self.create_relationship(aligned_pair_index=20)
        self.assertFalse(relationship.is_deletion())

        relationship = self.create_relationship(aligned_pair_index=2)
        self.assertFalse(relationship.is_deletion())

        relationship = self.create_relationship(aligned_pair_index=106)
        self.assertTrue(relationship.is_deletion())


    def test_is_mismatch(self):
        relationship = self.create_relationship(aligned_pair_index=20)
        self.assertFalse(relationship.is_mismatch())

        relationship = self.create_relationship(aligned_pair_index=3)
        self.assertTrue(relationship.is_mismatch())

        relationship = self.create_relationship(aligned_pair_index=106)
        self.assertTrue(relationship.is_mismatch())

        relationship = self.create_relationship(aligned_pair_index=0)
        self.assertTrue(relationship.is_mismatch())

        relationship = self.create_relationship(aligned_pair_index=1)
        self.assertTrue(relationship.is_mismatch())


    def test_is_between_pam_and_n20_for_ccn(self):

        relationship = self.create_relationship(aligned_pair_index=100, is_ngg=False)
        self.assertTrue(relationship.is_between_pam_and_n20())

        relationship = self.create_relationship(aligned_pair_index=86, is_ngg=False)
        self.assertFalse(relationship.is_between_pam_and_n20())

        relationship = self.create_relationship(aligned_pair_index=110, is_ngg=False)
        self.assertFalse(relationship.is_between_pam_and_n20())

    def test_is_between_pam_and_n20_for_ngg(self):

        relationship = self.create_relationship(aligned_pair_index=100)
        self.assertTrue(relationship.is_between_pam_and_n20())

        relationship = self.create_relationship(aligned_pair_index=86)
        self.assertFalse(relationship.is_between_pam_and_n20())

        relationship = self.create_relationship(aligned_pair_index=110)
        self.assertTrue(relationship.is_between_pam_and_n20())

        relationship = self.create_relationship(aligned_pair_index=111)
        self.assertFalse(relationship.is_between_pam_and_n20())

    def test_next_to_mismatch_or_indel(self):
        relationship = self.create_relationship(aligned_pair_index=100)
        self.assertFalse(relationship.next_to_mismatch_or_indel())

        relationship = self.create_relationship(aligned_pair_index=105)
        self.assertTrue(relationship.next_to_mismatch_or_indel())

        relationship = self.create_relationship(aligned_pair_index=106)
        self.assertFalse(relationship.next_to_mismatch_or_indel())

        relationship = self.create_relationship(aligned_pair_index=107)
        self.assertTrue(relationship.next_to_mismatch_or_indel())

        relationship = self.create_relationship(aligned_pair_index=5)
        self.assertTrue(relationship.next_to_mismatch_or_indel())

class TestDenotationIndex(unittest.TestCase):

    def test_create_denotation(self):
        denotation_A = DenotationIndex(index=5, cutsite=False)
        denotation_B = DenotationIndex(index=10, cutsite=True)

        self.assertEqual(denotation_A.index, 5)
        self.assertFalse(denotation_A.cutsite)
        self.assertEqual(denotation_A.representation, "|")

        self.assertEqual(denotation_B.index, 10)
        self.assertTrue(denotation_B.cutsite)
        self.assertEqual(denotation_B.representation, "||")

class TestCas9Denotations(unittest.TestCase):

    def test_apply_to_ccn_presentation_array(self):
        reference = ['A', 'T', 'T', 'A', '_', '_', 'G', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'T', 'T', 'T', 'C', 'A', 'G', 'C', 'A', 'A', 'G', 'G', 'C', '-', 'C', 'T', 'A', '-', '-', '-', '-', '-', '-', 'T', 'C', 'T', '_', 'G', '-', '-', '-', 'C', 'C', 'A', '-', '-', '-', '-', 'C', 'A', 'G', 'T', '-', '-', '-', '-', '-', 'G', 'A', 'G', '-', '-', 'G', 'T', 'C', '-', '-', '-', 'T', '_', '_', '_', '_', '_', '_', '_', '_', '_', 'G', 'A', 'C', 'T', 'C', 'C', 'C']

        read = ['G', 'T', 'T', 'C', 'G', 'C', 'G', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'C', 'C', 'T', 'A', 'A', 'G', 'G', '_', 'A', 'G', 'T', 'T', 'T', 'C', 'C', 'A', 'G', 'A', 'A', 'G', 'G', 'C', '-', 'C', 'C', 'A', '-', '-', '-', '-', '-', '-', 'T', 'T', 'T', 'G', 'G', '-', '-', '-', 'C', 'A', 'A', '-', '-', '-', '-', 'C', 'G', 'T', 'T', '-', '-', '-', '-', '-', 'G', 'G', 'G', '-', '-', 'G', 'C', 'C', '-', '-', '-', 'T', 'G', 'C', 'C', 'C', 'C', 'C', 'C', 'G', 'G', 'G', 'A', 'T', 'T', 'G', 'G', 'G']

        cutsite_index = 179
        pam_index = 173
        n20_pam_index = 176
        n20_index = 195
        aligned_pairs = [(0, 32), (1, 33), (2, 34), (3, 35), (4, None), (5, None), (6, 36), (7, 37), (8, 38), (9, 39), (10, 40), (11, 41), (12, 42), (13, 43), (14, 44), (15, 45), (16, 46), (17, 47), (18, 48), (19, 49), (20, 50), (21, 51), (22, 52), (23, 53), (24, 54), (25, 55), (26, 56), (27, 57), (28, 58), (29, 59), (30, 60), (31, 61), (32, 62), (33, 63), (34, 64), (35, 65), (36, 66), (37, 67), (38, 68), (39, 69), (40, 70), (41, 71), (42, 72), (43, 73), (44, 74), (45, 75), (46, 76), (47, 77), (48, 78), (49, 79), (50, 80), (51, 81), (52, 82), (53, 83), (54, 84), (55, 85), (56, 86), (57, 87), (58, 88), (59, 89), (60, 90), (61, 91), (62, 92), (63, 93), (64, 94), (65, 95), (66, 96), (67, 97), (68, 98), (69, 99), (70, 100), (71, 101), (72, 102), (73, 103), (74, 104), (75, 105), (76, 106), (77, 107), (78, 108), (79, 109), (80, 110), (81, 111), (82, 112), (83, 113), (84, 114), (85, 115), (86, 116), (87, 117), (88, 118), (89, 119), (90, 120), (91, 121), (92, 122), (93, 123), (94, 124), (95, 125), (96, 126), (97, 127), (98, 128), (99, 129), (100, 130), (101, 131), (102, 132), (103, 133), (104, 134), (105, 135), (106, 136), (107, 137), (108, 138), (109, 139), (110, 140), (111, 141), (112, 142), (113, 143), (114, 144), (115, 145), (116, 146), (117, 147), (118, 148), (119, 149), (120, 150), (121, 151), (122, 152), (123, 153), (124, 154), (125, 155), (126, 156), (127, 157), (128, 158), (129, 159), (130, 160), (131, 161), (132, 162), (133, 163), (134, 164), (135, 165), (136, 166), (137, 167), (138, 168), (139, 169), (140, 170), (141, 171), (142, 172), (143, 173), (144, 174), (145, 175), (146, 176), (147, 177), (148, 178), (149, 179), (None, 180), (150, 181), (151, 182), (152, 183), (153, 184), (154, 185), (155, 186), (156, 187), (157, 188), (158, 189), (159, 190), (160, 191), (161, 192), (162, 193), (163, 194), (164, 195), (165, 196), (166, 197), (167, 198), (168, 199), (169, 200), (170, 201), (171, 202), (172, 203), (173, 204), (174, 205), (175, 206), (176, 207), (177, None), (178, 208), (179, 209), (180, 210), (181, 211), (182, 212), (183, 213), (184, 214), (185, 215), (186, 216), (187, 217), (188, 218), (189, 219), (190, 220), (191, 221), (192, 222), (193, 223), (194, 224), (195, 225), (196, 226), (197, 227), (198, 228), (199, 229), (200, 230), (201, 231), (202, 232), (203, 233), (204, 234), (205, 235), (206, 236), (207, 237), (208, 238), (209, 239), (210, None), (211, None), (212, None), (213, None), (214, None), (215, None), (216, None), (217, None), (218, None), (219, 240), (220, 241), (221, 242), (222, 243), (223, 244), (224, 245), (225, 246)]
        is_ngg = False

        denotations = Cas9Denotations(cutsite_index, pam_index, n20_pam_index, n20_index, aligned_pairs, is_ngg)

        reference_string, read_string = denotations.apply_to_presentation(reference, read)

        expected_reference = "ATTA__G----------------------------------------------------------------------------------------------------------------------------------------|CCT|AAG||GAAGTTTCAGCAAGGC|-CTA------TCT_G---CCA----CAGT-----GAG--GTC---T_________GACTCCC"
        expected_read = "GTTCGCG----------------------------------------------------------------------------------------------------------------------------------------|CCT|AAG||G_AGTTTCCAGAAGGC|-CCA------TTTGG---CAA----CGTT-----GGG--GCC---TGCCCCCCGGGATTGGG"

        self.assertEqual(reference_string, expected_reference)
        self.assertEqual(read_string, expected_read)

    def test_apply_to_ngg_presentation_array(self):
        reference = ['G', 'A', 'C', '_', '_', '_', 'G', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'T', 'A', 'G', 'A', 'G', 'G', 'G', 'G', 'C', 'C', 'G', 'A', 'C', 'G', 'G', 'A', 'G', 'A', 'T', 'T', 'A', 'G', 'G', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']
        read = ['T', 'C', 'C', 'G', 'C', 'A', 'G', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'T', 'A', 'G', 'A', 'G', 'G', 'G', 'G', 'C', 'C', 'G', 'A', 'C', 'G', 'G', 'A', 'G', 'A', '_', 'T', 'A', 'G', 'G', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-']

        cutsite_index = 115
        pam_index = 121
        n20_pam_index = 118
        n20_index = 98

        aligned_pairs = [(0, 14), (1, 15), (2, 16), (3, None), (4, None), (5, None), (6, 17), (7, 18), (8, 19), (9, 20), (10, 21), (11, 22), (12, 23), (13, 24), (14, 25), (15, 26), (16, 27), (17, 28), (18, 29), (19, 30), (20, 31), (21, 32), (22, 33), (23, 34), (24, 35), (25, 36), (26, 37), (27, 38), (28, 39), (29, 40), (30, 41), (31, 42), (32, 43), (33, 44), (34, 45), (35, 46), (36, 47), (37, 48), (38, 49), (39, 50), (40, 51), (41, 52), (42, 53), (43, 54), (44, 55), (45, 56), (46, 57), (47, 58), (48, 59), (49, 60), (50, 61), (51, 62), (52, 63), (53, 64), (54, 65), (55, 66), (56, 67), (57, 68), (58, 69), (59, 70), (60, 71), (61, 72), (62, 73), (63, 74), (64, 75), (65, 76), (66, 77), (67, 78), (68, 79), (69, 80), (70, 81), (71, 82), (72, 83), (73, 84), (74, 85), (75, 86), (76, 87), (77, 88), (78, 89), (79, 90), (80, 91), (81, 92), (82, 93), (83, 94), (84, 95), (85, 96), (86, 97), (87, 98), (88, 99), (89, 100), (90, 101), (91, 102), (92, 103), (93, 104), (94, 105), (95, 106), (96, 107), (97, 108), (98, 109), (99, 110), (100, 111), (101, 112), (102, 113), (103, 114), (104, 115), (105, 116), (None, 117), (106, 118), (107, 119), (108, 120), (109, 121), (110, 122), (111, 123), (112, 124), (113, 125), (114, 126), (115, 127), (116, 128), (117, 129), (118, 130), (119, 131), (120, 132), (121, 133), (122, 134), (123, 135), (124, 136), (125, 137), (126, 138), (127, 139), (128, 140), (129, 141), (130, 142), (131, 143), (132, 144), (133, 145), (134, 146), (135, 147), (136, 148), (137, 149), (138, 150), (139, 151), (140, 152), (141, 153), (142, 154), (143, 155), (144, 156), (145, 157), (146, 158), (147, 159), (148, 160), (149, 161), (150, 162), (151, 163), (152, 164), (153, 165), (154, 166), (155, 167), (156, 168), (157, 169), (158, 170), (159, 171), (160, 172), (161, 173), (162, 174), (163, 175), (164, 176), (165, 177), (166, 178), (167, 179), (168, 180), (169, 181), (170, 182), (171, 183), (172, 184), (173, 185), (174, 186), (175, 187), (176, 188), (177, 189), (178, 190), (179, 191), (180, 192), (181, 193), (182, 194), (183, 195), (184, 196), (185, 197), (186, 198), (187, 199), (188, 200), (189, 201), (190, 202), (191, 203), (192, 204), (193, 205), (194, 206), (195, 207), (196, 208), (197, 209), (198, 210), (199, 211), (200, 212), (201, 213), (202, 214), (203, 215), (204, 216), (205, 217), (206, 218), (207, 219), (208, 220), (209, 221), (210, 222), (211, 223), (212, 224), (213, 225), (214, 226), (215, 227), (216, 228), (217, 229), (218, 230), (219, 231), (220, 232), (221, 233), (222, 234), (223, 235), (224, 236), (225, 237)]
        is_ngg = True

        denotations = Cas9Denotations(cutsite_index, pam_index, n20_pam_index, n20_index, aligned_pairs, is_ngg)

        reference_string, read_string = denotations.apply_to_presentation(reference, read)

        expected_reference = "GAC___G---------------------------------------------------------------------------------|TAGAGGGGCCGACGGAG||ATT|AGG|--------------------------------------------------------------------------------------------------------------------"
        expected_read = "TCCGCAG---------------------------------------------------------------------------------|TAGAGGGGCCGACGGAG||A_T|AGG|--------------------------------------------------------------------------------------------------------------------"

        self.assertEqual(reference_string, expected_reference)
        self.assertEqual(read_string, expected_read)

class TestPresenter(unittest.TestCase):

    def test_foo(self):
        pass

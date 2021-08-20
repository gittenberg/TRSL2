import unittest
from translation import MRNA, parameters


class TestMRNAMethods(unittest.TestCase):
    def test_attach_ribosome_at_start(self):
        testMRNA = MRNA.MRNA(index=0, length=parameters.mRNA_av_length, geneID=None, ribosomes={})
        testMRNA.attach_ribosome_at_start()
        self.assertEqual(testMRNA.ribosomes, {0: None}, msg="test_attach_ribosome_at_start...")

    def test_detach_ribosome(self):
        testMRNA = MRNA.MRNA(index=1, length=parameters.mRNA_av_length, geneID=None, ribosomes={})
        testMRNA.attach_ribosome_at_start()
        testMRNA.detach_ribosome(pos=0)
        self.assertEqual(testMRNA.ribosomes, {})

    def test_first_position_occupied(self):
        testMRNA = MRNA.MRNA(index=2, length=300, geneID=None, ribosomes={3: None})
        self.assertTrue(expr=testMRNA.first_position_occupied())
        testMRNA.detach_ribosome(pos=3)
        self.assertFalse(expr=testMRNA.first_position_occupied())

    def test_next_range_free(self):
        testMRNA = MRNA.MRNA(index=3, length=600, geneID=None, ribosomes={99: None, 300: None})
        self.assertTrue(expr=testMRNA.next_range_free(pos=12, by=24))
        self.assertFalse(expr=testMRNA.next_range_free(pos=90, by=300))

    def test_find_max_free_range(self):
        testMRNA = MRNA.MRNA(index=4, length=600, geneID=None, ribosomes={99: None, 300: None})
        self.assertEqual(first=testMRNA.find_max_free_range(pos=420), second=180)
        self.assertEqual(first=testMRNA.find_max_free_range(pos=297), second=3)

    def test_termination_condition(self):
        testMRNA1 = MRNA.MRNA(index=5, length=300, geneID=None, ribosomes={99: None, 300: None})
        self.assertTrue(expr=testMRNA1.termination_condition())
        testMRNA2 = MRNA.MRNA(index=6, length=330, geneID=None, ribosomes={99: None, 300: None})
        self.assertFalse(expr=testMRNA2.termination_condition())

    def test_prev_range_free(self):
        testMRNA = MRNA.MRNA(index=7, length=600, geneID=None, ribosomes={100: None, 300: None})
        self.assertTrue(expr=testMRNA.prev_range_free(pos=300, by=20))
        self.assertFalse(expr=testMRNA.prev_range_free(pos=300, by=220))

    def test_translocate_ribosome(self):
        testMRNA = MRNA.MRNA(index=7, length=600, geneID=None, ribosomes={100: None, 300: None})
        testMRNA.translocate_ribosome(pos=100, by=6)
        self.assertEqual(first=testMRNA.ribosomes, second={106: None, 300: None})

if __name__ == '__main__':
    unittest.main(verbosity=2)

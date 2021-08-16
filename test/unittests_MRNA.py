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
        testMRNA = MRNA.MRNA(index=3, length=600, geneID=None, ribosomes={99: None, 300: None})
        self.assertEqual(first=testMRNA.find_max_free_range(pos=420), second=180)
        self.assertEqual(first=testMRNA.find_max_free_range(pos=297), second=3)


if __name__ == '__main__':
    unittest.main(verbosity=2)
    #test
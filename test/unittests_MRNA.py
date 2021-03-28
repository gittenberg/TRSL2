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


if __name__ == '__main__':
    unittest.main(verbosity=2)

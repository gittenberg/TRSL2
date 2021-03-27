import unittest
from translation import MRNA, parameters


class TestMRNAMethods(unittest.TestCase):
    def test_attach_ribosome_at_start(self):
        testMRNA = MRNA.MRNA(index=0, length=parameters.mRNA_av_length, geneID=None, ribosomes={})
        testMRNA.attach_ribosome_at_start()
        self.assertEqual(testMRNA.ribosomes, {0: None})

    def test_detach_ribosome(self):
        testMRNA = MRNA.MRNA(index=0, length=parameters.mRNA_av_length, geneID=None, ribosomes={})
        testMRNA.attach_ribosome_at_start()
        testMRNA.detach_ribosome(pos=0)
        self.assertEqual(testMRNA.ribosomes, {})


if __name__ == '__main__':
    unittest.main()

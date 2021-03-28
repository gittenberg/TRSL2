import unittest
from translation import MRNA_specific, parameters


class TestMRNA_specificMethods(unittest.TestCase):
    def test_termination_condition(self):
        testMRNA = MRNA_specific.mRNA_spec(index=0, sequence="agcaaguga",
                                           geneID=None, ribosomes={6: None}, init_rate=1e-3)
        self.assertTrue(testMRNA.termination_condition())


if __name__ == '__main__':
    unittest.main(verbosity=3)

import unittest
import collections as col
from translation import MRNA, TRSL, parameters

class TestTRSLMethods(unittest.TestCase):
    def test_insert_tRNA(self):
        trsl = TRSL.TRSL(nribo=2, proteome=col.Counter({}), types_tRNA=2, n_mRNA=2, detail=False)
        # create mRNAs with unoccupied ribosomes
        trsl.mRNAs = [MRNA.MRNA(index=gene, length=9, ribosomes={3: None}) for gene in [0, 1]]
        # now insert tRNA into one of the mRNAs
        success = trsl.insert_tRNA(mRNA=trsl.mRNAs[0], pos=3, tRNA_type=1)
        # did it work overall
        self.assertTrue(success)
        # is the tRNA in the right place
        self.assertTrue(trsl.mRNAs[0].ribosomes=={3: 1})


if __name__ == '__main__':
    unittest.main(verbosity=2)

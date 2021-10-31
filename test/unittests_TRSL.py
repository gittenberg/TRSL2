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
        self.assertEqual(trsl.mRNAs[0].ribosomes, {3: 1})

    def test_release_tRNA(self):
        trsl = TRSL.TRSL(nribo=2, proteome=col.Counter({}), types_tRNA=2, n_mRNA=2, detail=False)
        # create mRNAs with unoccupied ribosomes
        trsl.mRNAs = [MRNA.MRNA(index=gene, length=9, ribosomes={3: None}) for gene in [0, 1]]
        # now insert tRNA into one of the mRNAs
        _ = trsl.insert_tRNA(mRNA=trsl.mRNAs[0], pos=3, tRNA_type=1)
        # finally release this tRNA
        success = trsl.release_tRNA(mRNA=trsl.mRNAs[0], pos=3, tRNA_type=1)
        # did it work overall
        self.assertTrue(success)
        # has the tRNA been removed
        self.assertEqual(trsl.mRNAs[0].ribosomes, {3: None})

    # not tested because obsolete # TODO
    def test_elongate_while_possible(self):
        self.assertTrue(True)

    # not tested because stochastic outcome
    def test_fill_empty_ribosomes(self):
        self.assertTrue(True)

    def test_elongate_one_step(self):
        trsl = TRSL.TRSL(nribo=2, proteome=col.Counter({}), types_tRNA=2, n_mRNA=2, detail=False)
        # create mRNAs with unoccupied ribosomes
        trsl.mRNAs = [MRNA.MRNA(index=gene, length=66, ribosomes={3: None}) for gene in [0, 1]]
        # now insert tRNA into one of the mRNAs
        _ = trsl.insert_tRNA(mRNA=trsl.mRNAs[0], pos=3, tRNA_type=1)
        success1 = trsl.elongate_one_step(trsl.mRNAs[0], current_pos=3)
        self.assertTrue(success1)
        trsl.mRNAs = [MRNA.MRNA(index=gene, length=12, ribosomes={3: None}) for gene in [0, 1]]
        # now insert tRNA into one of the mRNAs
        _ = trsl.insert_tRNA(mRNA=trsl.mRNAs[0], pos=3, tRNA_type=1)
        success2 = trsl.elongate_one_step(trsl.mRNAs[0], current_pos=3)
        # False because mRNA is too short
        self.assertFalse(success2)

    # TODO: unfinished!


if __name__ == '__main__':
    unittest.main(verbosity=2)

from translation import MRNA, parameters

testMRNA = MRNA.MRNA(index=0, length=parameters.mRNA_av_length, geneID=None, ribosomes={})
print(testMRNA.ribosomes)
testMRNA.attach_ribosome_at_start()
print(testMRNA.ribosomes)

# test
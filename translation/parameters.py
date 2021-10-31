#############################################################################################################
# used in MRNA_specific.py
#############################################################################################################
init_rate = 8.2e-07                        # average initiation rate in s^{-1}
stopcodons = ['uga', 'uaa', 'uag']
cr = 10                                    # ribosome footprint in codons

#############################################################################################################
# used in TRSL.py
#############################################################################################################
mRNA_av_length = 1251                      # source: http://bionumbers.hms.harvard.edu/bionumber.aspx?id=107678
types_tRNA = 42                            # number of types of tRNAs, including termination factor
n_mRNA = 60000                             # 60000 # number of mRNAs
                                           # (or 20000: http://book.bionumbers.org/how-many-mrnas-are-in-a-cell/)
nribo = 200000                             # number of ribosomes
                                           # http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=100267&ver=13&trm=ribosomes/cell
init_rate_low = 3.5e-06                    # initiation probability at mRNA 5' end (lower bound)
init_rate_high = 0.115                     # initiation probability at mRNA 5' end (upper bound)
V = 4.2e-17                                # m^3 # cell volume
                                           # http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000865
n_genes = 3795                             # number of protein-coding genes in the experiment
n_tRNA = 3300000                           # number of tRNAs
                                           # http://nar.oxfordjournals.org/content/suppl/2011/04/23/gkr300.DC1/Supplemental_File_S2.pdf gives 2984788
avogadro = 6.022e23
lambda_tRNA = 1.5e-8                       # m # characteristic length tRNA
D_tRNA = 8.42e-11                          # m^2/s # diffusion coefficient of tRNA
tau_tRNA = lambda_tRNA ** 2 / 6. / D_tRNA  # s # char. time for tRNA
num_pos_tRNA = V / lambda_tRNA ** 3        # number of discrete positions for tRNA
lambda_ribo = 3e-8                         # m # characteristic length ribosomes
D_ribo = 3e-13                             # m^2/s # diffusion coefficient of ribosomes
tau_ribo = lambda_ribo ** 2 / 6. / D_ribo  # s # char. time for ribosomes
num_pos_ribo = V / lambda_ribo ** 3        # number of discrete positions for ribosomes
competition = 7.78e-4                      # tRNA competition coefficient

#############################################################################################################
# used in TRSL_specific.py
#############################################################################################################
'''
these codons are on the mRNA from 5' to 3'
non-primary source also: http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002203.s002&type=supplementary
'''
genetic_code = {
    'uuu': 'F', 'ucu': 'S', 'uau': 'Y', 'ugu': 'C',
    'uuc': 'F', 'ucc': 'S', 'uac': 'Y', 'ugc': 'C',
    'uua': 'L', 'uca': 'S', 'uaa': '*', 'uga': '*',  # '*'==stop
    'uug': 'L', 'ucg': 'S', 'uag': '*', 'ugg': 'W',
    'cuu': 'L', 'ccu': 'P', 'cau': 'H', 'cgu': 'R',
    'cuc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
    'cua': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
    'cug': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
    'auu': 'I', 'acu': 'T', 'aau': 'N', 'agu': 'S',
    'auc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
    'aua': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
    'aug': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
    'guu': 'V', 'gcu': 'A', 'gau': 'D', 'ggu': 'G',
    'guc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
    'gua': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
    'gug': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
}

'''
source: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC147333/
we use 5'-3' convention for codon and anticodon;
we ignore special nucleotides i and psi and use a and u instead;
some simplification by making this a 1:1 relationship
note that https://academic.oup.com/bioinformatics/article/28/18/i340/248836 are mapping acg to cau not cgu
'''

codon_anticodon = {
    'uuu': 'gaa', 'ucu': 'aga', 'uau': 'gua', 'ugu': 'gca',
    'uuc': 'gaa', 'ucc': 'aga', 'uac': 'gua', 'ugc': 'gca',
    'uua': 'uaa', 'uca': 'uga', 'uaa': '*', 'uga': '*',  # '*'==stop
    'uug': 'caa', 'ucg': 'cga', 'uag': '*', 'ugg': 'cca',
    'cuu': 'gag', 'ccu': 'agg', 'cau': 'gug', 'cgu': 'acg',
    'cuc': 'gag', 'ccc': 'agg', 'cac': 'gug', 'cgc': 'acg',
    'cua': 'uag', 'cca': 'ugg', 'caa': 'uug', 'cga': 'acg',
    'cug': 'uag', 'ccg': 'ugg', 'cag': 'cug', 'cgg': 'ccg',
    'auu': 'aau', 'acu': 'agu', 'aau': 'guu', 'agu': 'gcu',
    'auc': 'aau', 'acc': 'agu', 'aac': 'guu', 'agc': 'gcu',
    'aua': 'uau', 'aca': 'ugu', 'aaa': 'uuu', 'aga': 'ucu',
    'aug': 'cau', 'acg': 'cgu', 'aag': 'cuu', 'agg': 'ccu',
    'guu': 'aac', 'gcu': 'agc', 'gau': 'guc', 'ggu': 'gcc',
    'guc': 'aac', 'gcc': 'agc', 'gac': 'guc', 'ggc': 'gcc',
    'gua': 'uac', 'gca': 'ugc', 'gaa': 'uuc', 'gga': 'ucc',
    'gug': 'cac', 'gcg': 'ugc', 'gag': 'cuc', 'ggg': 'ccc'
}

''' source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159466/bin/supp_39_15_6705__index.html
    https://academic.oup.com/nar/article/39/15/6705/1022014/The-role-of-tRNA-and-ribosome-competition-in
    anticodons are in 5'-3' direction (same convention as above)'''
tRNA_types = {
    1: {'anticodon': 'ugc', 'abundancy': 55351},  # reverse complement the anticodon to look it up
    2: {'anticodon': 'agc', 'abundancy': 121771},
    3: {'anticodon': 'ucu', 'abundancy': 121771},
    4: {'anticodon': 'ccu', 'abundancy': 11070},
    5: {'anticodon': 'ccg', 'abundancy': 11070},
    6: {'anticodon': 'acg', 'abundancy': 66421},
    7: {'anticodon': 'guu', 'abundancy': 110701},
    8: {'anticodon': 'guc', 'abundancy': 177122},
    9: {'anticodon': 'gca', 'abundancy': 44280},
    10: {'anticodon': 'uug', 'abundancy': 88561},
    11: {'anticodon': 'cug', 'abundancy': 11070},
    12: {'anticodon': 'uuc', 'abundancy': 154982},
    13: {'anticodon': 'ucc', 'abundancy': 33210},
    14: {'anticodon': 'ccc', 'abundancy': 22140},
    15: {'anticodon': 'gcc', 'abundancy': 177122},
    16: {'anticodon': 'gug', 'abundancy': 77491},
    17: {'anticodon': 'uau', 'abundancy': 22140},
    18: {'anticodon': 'aau', 'abundancy': 143911},
    19: {'anticodon': 'uag', 'abundancy': 33210},
    20: {'anticodon': 'gag', 'abundancy': 11070},
    21: {'anticodon': 'uaa', 'abundancy': 77491},
    22: {'anticodon': 'caa', 'abundancy': 110701},
    23: {'anticodon': 'uuu', 'abundancy': 77491},
    24: {'anticodon': 'cuu', 'abundancy': 154982},
    25: {'anticodon': 'cau', 'abundancy': 55351},  # 26 does not seem to exist
    27: {'anticodon': 'gaa', 'abundancy': 110701},
    28: {'anticodon': 'agg', 'abundancy': 22140},
    29: {'anticodon': 'ugg', 'abundancy': 110701},
    30: {'anticodon': 'gcu', 'abundancy': 33210},
    31: {'anticodon': 'uga', 'abundancy': 33210},
    32: {'anticodon': 'aga', 'abundancy': 121771},
    33: {'anticodon': 'cga', 'abundancy': 11070},
    34: {'anticodon': 'ugu', 'abundancy': 44280},
    35: {'anticodon': 'agu', 'abundancy': 121771},
    36: {'anticodon': 'cgu', 'abundancy': 11070},
    37: {'anticodon': 'cca', 'abundancy': 66421},
    38: {'anticodon': 'gua', 'abundancy': 88561},
    39: {'anticodon': 'uac', 'abundancy': 22140},
    40: {'anticodon': 'aac', 'abundancy': 154982},
    41: {'anticodon': 'cac', 'abundancy': 22140},
    42: {'anticodon': '*', 'abundancy': 18000},  # termination factor
    43: {'anticodon': 'cuc', 'abundancy': 22140}
}

anticodon_index = {tRNA_types[i]['anticodon']: i for i in tRNA_types}

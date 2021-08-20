cr = 10                                    # ribosome footprint in codons
mRNA_av_length = 1251                      # source: http://bionumbers.hms.harvard.edu/bionumber.aspx?id=107678

# used in MRNA_specific.py
init_rate = 8.2e-07                        # average initiation rate in s^{-1}
stopcodons = ['uga', 'uaa', 'uag']

# used in TRSL.py
nribo = 200000                             # number of ribosomes
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

import collections as col
import pickle

import translation.MRNA_specific
import translation.TRSL_specific

print("TRSL_setup: starting")

examplesequence_1 = "ggg uuu uca uca uuu gag gac gau gua aaa ggg uaa".replace(' ', '')
examplesequence_2 = "aug aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa cug ccc gag ggg uuu uca uca uuu gag gac aaa\
 cug ccc gag ggg uuu uca uca uaa".replace(' ', '')

conf = {1: {
    # mini test case run
    'exome': {1: examplesequence_1, 2: examplesequence_2},
    'transcriptome': {1: 2, 2: 1},
    'init_rates': {1: 0.1, 2: 0.1},
    'description': 'test configuration with 2 genes in 3 transcripts',
    'tRNA': col.Counter({i: 100 for i in translation.TRSL_specific.tRNA_types})},
    2: {
    'exome': pickle.load(open("./parameters/orf_coding.p", "rb")),
    'transcriptome': pickle.load(open("./parameters/transcriptome_plotkin_20000.p", "rb")),
    'init_rates': pickle.load(open("./parameters/init_rates_plotkin.p", "rb")),
    'description': '20000 transcriptome, full exome, no decay, Plotkin initiation probabilities'},
    8: {
    # base case run with reduced number of genes (limited by availability of initiation rates)
    'exome': pickle.load(open("./parameters/orf_coding.p", "rb")),
    'transcriptome': pickle.load(open("./parameters/transcriptome_shah.p", "rb")),
    'init_rates': pickle.load(open("./parameters/init_rates_plotkin.p", "rb")),
    'description': 'updated Shah transcriptome, full exome, no decay, updated initiation rates according to Shah'}
}

if __name__ == "__main__":
    for i in [8]:  # set configuration_id
        genes = list(set(conf[i]['exome']) & set(conf[i]['transcriptome']) & set(conf[i]['init_rates']))
        conf[i]['decay_constants'] = None
        print("found %s genes in common." % len(genes))

        mRNAs = []
        counter = 0
        for gene in genes:
            if conf[i]['init_rates']:
                if gene in conf[i]['transcriptome'] and gene in conf[i]['init_rates']:
                    # print "abundances and initiation rates available for gene:", gene
                    for instance in range(conf[i]['transcriptome'][gene]):
                        mRNAs.append(translation.MRNA_specific.mRNA_spec(index=counter,
                                                                         sequence=conf[i]['exome'][gene],
                                                                         geneID=gene,
                                                                         ribosomes={},
                                                                         init_rate=conf[i]['init_rates'][gene]))
                        # do not just multiply the list
                        counter += 1
            else:
                if gene in conf[i]['transcriptome']:
                    print("abundancies but no initiation rate available for gene:", gene)
                    for instance in range(conf[i]['transcriptome'][gene]):
                        mRNAs.append(translation.MRNA_specific.mRNA_spec(index=counter,
                                                                         sequence=conf[i]['exome'][gene],
                                                                         geneID=gene,
                                                                         ribosomes={}))  # do not just multiply the list
                        counter += 1
        print("built gene library, next: run TRSL_spec.")

        description = conf[i]['description']
        print(description)

        duration = 60.0

        tr = translation.TRSL_specific.TRSL_spec(mRNAs, conf[i]['exome'], conf[i]['decay_constants'],
                                                 nribo=20,
                                                 detail=True)
        # overwrite tRNA:
        if 'tRNA' not in conf[i]:
            tr._tRNA = col.Counter({i: translation.TRSL_specific.tRNA_types[i]['abundancy']
                                    for i in translation.TRSL_specific.tRNA_types})  # full tRNA
        else:
            tr._tRNA = conf[i]['tRNA']
        # print tr._tRNA
        tr._tRNA_free = col.Counter({i: int(tr._tRNA[i]) for i in translation.TRSL_specific.tRNA_types})
        # tRNA not bound to ribosomes
        tr._tRNA_bound = tr._tRNA - tr._tRNA_free
        # tRNA bound to ribosomes

        print("found {} mRNAs".format(len(mRNAs)))
        print("found {} genes".format(len(conf[i]['exome'])))
        print("found {} tRNA molecules".format(sum(tr._tRNA.values())))
        print("solving...")

        tr.solve_internal(0.0, duration, deltat=0.05)

        tr.dump_results(description, dirname="")

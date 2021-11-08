"""
protein translation as a probabilistic, discrete model (generic version)

This module generates an dynamical system of polysome translating proteins.
The proteins are generic (not sequence-specific).

Required steps and energy consumption:

1. initiation (cost: 2 GTP)
1.1 ribosome attachment (first subunit)
1.2 AA-tRNA binding at P site.
1.3 ribosome attachment (second subunit)
   initiation complex:  2 GTP -> 2 GDP
   We currently merge 1.1, 1.2, 1.3 into one binomial experiment

2. elongation
2.1 AA-tRNA binding at A site (cost: (n-1) GTP -> (n-1) GDP)
2.2 peptide bond formation and peptide elongation (no cost for bond)
2.3 translocation P>E, A>P (cost: (n-1) GTP -> (n-1) GDP)
2.4 release tRNA from E
    We carry out 2.4 (release) before 2.1 (binding) so we do not have to distinguish E, P, A sites
    Total: 2 (n-1) GTP -> 2 (n-1) GDP

3. termination (cost: 2 GTP -> 2 GDP, in reality, one of the GTPs is an ATP but we want to avoid using ADP in the model)

Total energy balance per protein:
We assume tRNA is activated outside of the model
2 (n+1) GTP -> 2 (n+1) GDP
"""

import sys
import cProfile
import pstats
import random as ran
import collections as col
import math
import logging as log
import os.path
import pickle

import numpy as np
import numpy.random as npr

import translation.MRNA

# read model parameters
from translation.parameters import *


class TRSL(object):
    """
    class representing a translational network
    """

    # initiation and auxiliary functions
    ###################################################################################################################
    def __init__(self, nribo=nribo, proteome=None, types_tRNA=types_tRNA, n_mRNA=n_mRNA,
                 mRNA_length=mRNA_av_length, detail=False):
        """
        initializes the parameters of the translation process
        """
        if proteome is None:
            proteome = col.Counter({})
        log.info("__init__: initializing TRSL")

        # Parameters
        ###############################################################################################################
        self.types_tRNA = types_tRNA

        # Initial values
        ##################################################################################################################################
        self.n_mRNA = n_mRNA
        self.ribo_free = nribo
        self.p_init = math.sqrt(init_rate_low * init_rate_high) # geometric mean btw the lowest and highest possible value

        self.GTP = 1e3 * avogadro * V  # GTP molecules (made up)
        self.GDP = 0                   # GDP molecules
        self.ATP = 1e3 * avogadro * V  # ATP molecules (made up)
        self.AMP = 0                   # AMP molecules

        self.timerange = []
        self.timecourses = {}
        self._tRNA = col.Counter({i: int(0.5 + n_tRNA / (self.types_tRNA * 1.0))
                                  for i in range(1, self.types_tRNA + 1)})
        # uniform distribution because translation in this file is not specific

        # Warning: if ribosomes are not passed explicitely in the following MRNA constructor, they will be passed by
        # reference and all MRNAs will share the same ribosome!
        self.mRNAs = [translation.MRNA.MRNA(index=gene, length=mRNA_length, ribosomes={})
                      for gene in [ran.randint(1, n_genes) for k in range(self.n_mRNA)]]  # randomized gene expressions

        self.proteins = proteome  # contains protein IDs and counts not including polypeptides in statu nascendi
        self.protein_length = sum(self.proteins.values())  # equals number of peptide bonds?

        self.init_rate = self.p_init / tau_ribo / num_pos_ribo   # 8.157e-07 s^-1 (yeast)
        self.elong_rate = competition / tau_tRNA / num_pos_tRNA  # 0.000140 s^-1  (yeast)

        self.detail = detail  # whether details are saved (e.g. ribosomes for every time step)

    @property
    def tRNA_bound(self):
        return self._tRNA_bound

    @tRNA_bound.setter
    def tRNA_bound(self, value):
        self._tRNA_bound = value

    @property
    def ribo_free(self):
        return self._ribo_free

    @ribo_free.setter
    def ribo_free(self, value):
        self._ribo_free = value

    @property
    def mRNAs(self):
        return self._mRNAs

    @mRNAs.setter
    def mRNAs(self, value):
        """
        This setter sets mRNA object and keeps the bound ribosomes and tRNA in sync.

        @:var value: list of mRNA objects
        """
        self.ribo_bound = sum([len(m.ribosomes) for m in value])  # all bound ribosomes
        all_ribos = self.ribo_bound + self.ribo_free
        # test = sum([len(m.ribosomes.keys()) for m in value])
        self.ribo_free = all_ribos - self.ribo_bound
        # free tRNAs for new set of polysomes
        trnas_in_polysomes = col.Counter()
        for m in value:
            m_bound_tRNA = m.ribosomes.values()
            for t in m_bound_tRNA:
                trnas_in_polysomes[t] += 1

        self.tRNA_free = self._tRNA - trnas_in_polysomes
        self.tRNA_bound = trnas_in_polysomes
        self._mRNAs = value

    @property
    def ribo_bound(self):
        return self._ribo_bound

    @ribo_bound.setter
    def ribo_bound(self, value):
        self._ribo_bound = value

    # def __getitem__(self, i): # TODO: loeschen?
    #     """
    #     This allows to address the modeldict like trsl['vars']
    #     """
    #     return self.modeldict[i]

    def inspect(self):
        """
        Print all dictionaries in TRSL
        """
        print("----------------------------------------------------------")
        for key in sorted(self.__dict__): # TODO: vereinfachen, wenn modeldict nicht mehr gebraucht
            if key != "modeldict":
                print(key, ":")
                print(self.__dict__[key])
                print("----------------------------------------------------------")
            else:
                print(key, ":")
                for subkey in self.__dict__[key]:
                    print("\t\t", subkey, ":")
                    # print "\t\t", self.__dict__[key][subkey]
                    print("----------------------------------------------------------")

    def get_state(self):
        """
        Get a dictionary of all the defining properties of the simulation.

        :return: dict
        """
        import time
        now = time.strftime("%Y%m%d_%H%M", time.gmtime())

        results = dict()
        results['proteome'] = self.proteins
        results['peptide_bonds'] = self.protein_length
        results['transcriptome'] = self.mRNAs
        results['timerange'] = self.timerange
        results["timecourses"] = self.timecourses
        results["description"] = ""
        results["time_stamp"] = now
        results["n_ribosomes"] = self.ribo_bound + self.ribo_free
        # results["collisions"] = self.collision      # TODO: needs to go to TRSL_specific
        # results["nocollisions"] = self.nocollision  # TODO: needs to go to TRSL_specific
        results["n_tRNA"] = sum(self._tRNA.values())
        duration = self.timerange[-1] - self.timerange[0] + 1
        results["duration"] = duration
        return results

    def dump_results(self, description='results', dirname=r'../results'):
        """
        Save results of the simulation to a pickle file in the ../results directory.
        The name is generated using the given description and a timestamp.

        @param description: readable string describing the simulation
        @param dirname: directory to save the results
        """
        results = self.get_state()
        results["description"] = description

        filename = "{}_{}_{}_ribosomes_{}s.p".format(description, results['time_stamp'], results["n_ribosomes"],
                                                     str(int(results["duration"])).zfill(4))
        pathname = os.path.join(dirname, filename)
        pickle.dump(results, open(pathname, "wb"))
        print(description)

    # functions used in simulation
    ####################################################################################################################

    def insert_tRNA(self, mRNA, pos, tRNA_type):
        """
        inserts a tRNA of type tRNA_type at position pos
        returns True iff successful
        """
        if pos in mRNA.ribosomes and self.tRNA_free[tRNA_type] >= 1 and not mRNA.ribosomes[pos]:
            # there has to be a ribosome at pos and there has to be tRNA of that type available and there cannot be a tRNA yet at pos on the position
            # log.debug("insert_tRNA: inserting tRNA %s on mRNA %s at position %s", tRNA_type, self.mRNAs.index(mRNA), pos)
            self.tRNA_free[tRNA_type] -= 1
            self.tRNA_bound[tRNA_type] += 1
            mRNA.ribosomes[pos] = tRNA_type  # tRNA now bound
            # log.debug("insert_tRNA: ribosomes: tRNA on mRNA %s are now %s", self.mRNAs.index(mRNA), mRNA.ribosomes)
            success = True
        elif not pos in mRNA.ribosomes:
            log.warning("insert_tRNA: cannot insert tRNA at %s: no ribosome present", pos)
            success = False
        elif mRNA.ribosomes[pos] != None:
            log.warning("insert_tRNA: cannot insert tRNA at %s: tRNA already bound", pos)
            success = False
        elif self.tRNA_free[tRNA_type] < 1:
            # log.warning("insert_tRNA: cannot insert tRNA type %s at %s: not enough free tRNA", tRNA_type, pos)
            # log.warning("insert_tRNA: tRNA_type is: %s", tRNA_type)
            # log.warning("insert_tRNA: bound tRNA is: %s", self.tRNA_bound)
            success = False
        else:
            log.warning("insert_tRNA: failed for unspecified reason")
            success = False
        return success

    def release_tRNA(self, mRNA, pos, tRNA_type):
        """
        releases a tRNA molecule of type tRNA_type from mRNA
        """
        if mRNA.ribosomes[pos] == tRNA_type and self._tRNA_bound[tRNA_type] > 0:
            self.tRNA_bound[tRNA_type] -= 1
            self.tRNA_free[tRNA_type] += 1
            mRNA.ribosomes[pos] = None
            #log.debug("release_tRNA: successfully released tRNA %s from pos %s", tRNA_type, pos)
            success = True
        else:
            #log.warning("release_tRNA: failed to release tRNA from pos %s", pos)
            #log.warning("release_tRNA: mRNA.ribosomes[pos] == %s", mRNA.ribosomes[pos])
            #log.warning("release_tRNA: self.tRNA_bound[tRNA_type] == %s", self._tRNA_bound[tRNA_type])
            success = False
        return success

    # TODO: Loeschkandidat
    # def elongate_while_possible(self, mRNA, k, current_pos):
    #     """
    #     attempts to elongate the protein on mRNA by at most k AAs starting at current_pos
    #     stops if steric hindrance by another ribosome , or end of mRNA is encountered
    #     NOTE: function is now obsolete
    #     """
    #     free_range = mRNA.find_max_free_range(current_pos)
    #     # log.debug("elongate_while_possible: found free range of %s nts downstream of %s", free_range, current_pos)
    #     codons = min(k, free_range / 3)  # number k of available tRNAs and sterically free codons limit elongation # integer division
    #     codons = max(0, codons)  # not negative
    #     codons = min((mRNA.length - current_pos) / 3, codons)  # cannot translate behind end of mRNA
    #     # log.debug("elongate_while_possible: %s tRNAs, free range of %s nts, trying to elongate by %s codons", k, free_range, codons)
    #     if self.GTP >= codons:
    #         if codons > 0:
    #             # log.debug("elongate_while_possible: possible to translocate by %s codons", codons)
    #             # elongation: release tRNA
    #             previous_type = mRNA.ribosomes[current_pos]  # type to be released at ribo_pos
    #             self.release_tRNA(mRNA, current_pos, previous_type)
    #             # translocation: move ribosome
    #             mRNA.translocate_ribosome(current_pos, by=3 * codons)
    #             # bind AA-tRNA
    #             last_type = ran.choice(self._tRNA.keys())  # type to be inserted at current_pos
    #             # log.debug("elongate_while_possible: last position was %s, attempting to insert tRNA at position %s", current_pos, current_pos+3*codons)
    #             self.insert_tRNA(mRNA, current_pos + 3 * codons, last_type)  # try to insert AA-tRNA in the ribosome
    #             # translocation: elongate proteinlength
    #             self.protein_length += codons
    #             self.GTP -= codons
    #             self.GDP += codons
    #     else:
    #         log.warning("elongate_while_possible: not enough GTP")
    #     # log.debug("elongate_while_possible: ribosomes: tRNA is now %s", mRNA.ribosomes)
    #     # log.debug("elongate_while_possible: protein length is now %s", self.proteinlength)

    def fill_empty_ribosomes(self, this_mRNA, time):
        """
        iterate through every empty ribosome and try to diffuse the required tRNA into the site
        returns True if at least one attempt was successful, else false
        """
        change_occurred = False
        # empty ribosomes on this particular mRNA
        # the ribosome at position 0 is filled during the initiation process (not modelled)
        empty_ribos = [key for key in this_mRNA.ribosomes if this_mRNA.ribosomes[key] is None and key != 0]
        # log.debug("fill_empty_ribosomes: this_mRNA.index = {}".format(this_mRNA.index))
        # log.debug("fill_empty_ribosomes: empty_ribos = {}".format(empty_ribos))
        # log.debug("fill_empty_ribosomes: this_mRNA.ribosomes = {}".format(this_mRNA.ribosomes))
        for ribo_pos in empty_ribos:
            required_tRNA_type = ran.choice(self._tRNA.keys())  # random type to be inserted
            tRNA_diffusion_probability = self.elong_rate * time  # ignoring wobble in the unspecific model
            failure_probability = (1 - tRNA_diffusion_probability) ** self.tRNA_free[required_tRNA_type]
            randomnumber = ran.random()  # TODO: try Poisson approximation if faster
            # log.debug("fill_empty_ribosomes: from mRNA %s: failure probability is %s at mRNA position %s. Available time: %s", this_mRNA, failure_probability, ribo_pos, time)
            success = not (randomnumber < failure_probability)  # this means the required tRNA type has diffused to the site
            # log.debug("fill_empty_ribosomes: success: {}".format(success))
            if success:
                # log.debug('fill_empty_ribosomes: matching tRNA diffused to initiation site')
                if not self.insert_tRNA(this_mRNA, ribo_pos, required_tRNA_type):
                    log.warning("fill_empty_ribosomes: unsuccessful attempt to insert tRNA")
                else:
                    # log.debug("fill_empty_ribosomes: successful attempt to insert tRNA")
                    change_occurred = True
        return change_occurred

    def elongate_one_step(self, mRNA, current_pos):
        """
        attempts to elongate the protein on mRNA by one AA at current_pos
        stops if steric hindrance by another ribosome or end of mRNA is encountered
        """
        # free_codons = (mRNA.find_max_free_range(current_pos) - 3 * translation.MRNA.cr) / 3  # integer division on purpose # TODO: old version - wrong?
        free_codons = mRNA.find_max_free_range(current_pos) / 3  # integer division on purpose
        if self.GTP >= 1 and free_codons > 0:
            # log.debug("elongate_one_step: possible to translocate")
            previous_type = mRNA.ribosomes[current_pos]  # type to be released at ribo_pos
            # log.debug("elongate_one_step: mRNA.ribosomes = %s", mRNA.ribosomes)
            # log.debug("elongate_one_step: self.tRNA_bound = %s", self.tRNA_bound)
            # translocation: move ribosome
            mRNA.translocate_ribosome(current_pos, by=3)
            self.release_tRNA(mRNA, current_pos+3, previous_type)
            # translocation: elongate proteinlength
            self.protein_length += 1
            self.GTP -= 2
            self.GDP += 2
            return True
        else:
            if free_codons <= 0:
                # log.warning("elongate_one_step: not possible: no free codon")
                return False
            else:
                # log.warning("elongate_one_step: not possible: not enough GTP or other reason")
                return False
            # log.debug("elongate_one_step: ribosomes: tRNA is now %s", mRNA.ribosomes)
            # log.debug("elongate_one_step: protein length is now %s", self.proteinlength)

    def elongate_mRNA(self, mRNA):
        """
        translocates all ribosomes on mRNA by one step where possible
        """
        # log.debug("elongate_mRNA: from mRNA %s: ribosomes on this mRNA are: %s", mRNA.index, mRNA.ribosomes)
        occupied_ribos = [key for key in mRNA.ribosomes if (mRNA.ribosomes[key] is not None or key==0)] # the first codon is always occupied by tRNA^Met_i
        for ribo_pos in occupied_ribos:  # TODO: test reverse list and other sort orders
            self.elongate_one_step(mRNA, ribo_pos)

    # functions used for process control
    ##################################################################################################################################

    def update_initiation(self, deltat, mRNA):
        """
        performs binomial random experiment to attach ribosome at start codon

        :param deltat: duration parameter driving the initiation probability
        :param mRNA:   mRNA object to which ribosome is attached
        :return:
        """
        # log.info("update_initiation: starting")
        # log.debug("update_initiation: found mRNA %s", j)
        k = npr.binomial(self.ribo_free, self.init_rate * deltat, 1)[0]  # number of ribosomes that diffuse to the initiation site during deltat
        # log.debug("update_initiation: %s ribosomes diffused to init site at mRNA %s", k, mRNA.index)
        for i in range(k):  # currently k>1 will not attach k ribosomes, TODO:
            if not mRNA.first_position_occupied():
                # log.debug("update_initiation: found mRNA with free first position")
                if self.GTP > 0 and self.ATP > 0:
                    if mRNA.attach_ribosome_at_start():
                        log.debug("update_initiation: attaching ribosome at start of mRNA %s", mRNA.index)
                        self.ribo_bound += 1
                        self.ribo_free -= 1
                        self.GTP -= 2
                        self.GDP += 2
                    else:
                        log.warning("update_initiation: unsuccessful attempt to attach ribosome")
                else:
                    log.warning("update_initiation: no GTP or no ATP")
            else:
                # log.warning("update_initiation: unsuccessful attempt to attach ribosome: first position occupied")
                pass

    def update_elongation(self, deltat, mRNA):
        """
        while a change still occurs:
            fill all empty ribosomes by tRNA diffusion
            if possible:
                elongate all polypeptides
                all occupied ribosomes move by one step and lose bound tRNA
            halve time interval and continue
        """
        log.info("update_elongation: starting mRNA %s, geneID %s", mRNA.index, mRNA.geneID)
        change_flag = True
        available_time = deltat
        while change_flag:  # while there is a change in tRNA or ribosome position
            change_flag = self.fill_empty_ribosomes(mRNA, available_time)  # if a tRNA bound this becomes True
            self.elongate_mRNA(mRNA)  # translocate all ribosomes by one step if possible
            available_time *= 0.5
            log.debug("update_elongation: from mRNA %s: halving time, available time is now %s", mRNA.index, available_time)

    def update_termination(self, mRNA, time):
        """
        check for termination condition
        if fulfilled, detach ribosome and increase number of proteins
        """
        #log.info("update_termination: starting")
        if self.GTP >= 1:
            if mRNA.termination_condition():
                # log.debug("update_termination: mRNA.ribosomes = %s", mRNA.ribosomes)
                release_pos = max(mRNA.ribosomes.keys())
                # release_type = mRNA.ribosomes[release_pos]
                # log.debug("update_termination: terminating translation at mRNA %s, release position %s, release type %s", mRNA.index, release_pos, release_type)
                # self.release_tRNA(mRNA, release_pos, release_type)  # commented because there is no tRNA at the stop codon
                mRNA.detach_ribosome(release_pos)
                self.ribo_bound -= 1
                self.ribo_free += 1
                if not mRNA.geneID in self.proteins:
                    self.proteins[mRNA.geneID] = 1  # add first protein of type mRNA.geneID
                else:
                    self.proteins[mRNA.geneID] += 1  # add another protein of type mRNA.geneID
                self.GTP -= 2  # in the cell one of these is an ATP but this is equivalent
                self.GDP += 2
                if mRNA.tic and mRNA.toc:
                    mRNA.toc -= 1                             # one ribosome falls off the mRNA
                if mRNA.tic and mRNA.toc==1:                  # time measurement ongoing and the last ribosome just fell off
                    mRNA.tic_toc.append((mRNA.tic, time))
                    mRNA.tic = False                          # reset time measurement
                    mRNA.toc = False
        else:
            log.warning("update_termination: not enough GTP")
    
    def update_processes(self, deltat, time):
        """
        for all mRNAs update initiation, elongation and termination
        """
        for mRNA in self.mRNAs:
            if mRNA.ribosomes != {}:  # mRNAs without ribosomes do not need the following
                self.update_termination(mRNA, time)
                self.update_elongation(deltat, mRNA)
            self.update_initiation(deltat, mRNA)
        # self.update_protein_decay(deltat)
        
    def initialize_solve_internal(self, fieldnames):
        self.timecourses = {fieldname: [] for fieldname in fieldnames}

        # tRNAs are separate because they are many:
        for tRNA_type in self.tRNA_free:
            self.timecourses["tRNA_free_" + str(tRNA_type).zfill(2)] = []

    def update_solve_internal(self, deltat, fields, start, time):
        # update standard
        for fieldname, fieldvalue in fields.items():
            self.timecourses[fieldname].append(fieldvalue)

        # now update proteins
        for gene_id in self.proteins:
            if gene_id in self.timecourses:  # if there is already protein of this ID
                # then append the protein count
                self.timecourses[gene_id].append(self.proteins[gene_id])
            else:  # if this is the first time protein of this ID terminated
                # create zeros for the first (time-start)/deltat values
                self.timecourses[gene_id] = [0] * int((time - start) / deltat)
                # only then append the protein count
                self.timecourses[gene_id].append(self.proteins[gene_id])

        # now update specific tRNA_free
        for tRNA_type in self.tRNA_free:
            self.timecourses["tRNA_free_" + str(tRNA_type).zfill(2)].append(self.tRNA_free[tRNA_type])
            # log.info("solve_internal: tRNA_free type %s: %s molecules", tRNA_type, self.tRNA_free[tRNA_type])

    def write_last_polysome(self, description):
        # this is a legacy of the old version when all polysomes were saved (now only last polysome)
        # if detail option, then initiate timecourse for every polysome
        if self.detail:
            import shelve
            import time
            import copy

            now = time.strftime("%Y%m%d_%H%M", time.gmtime())
            timestamp = now
            ribosomes_database = shelve.open('../results/ribosome_timecourses_{}_{}.db'.format(description, timestamp), writeback=True)
            for mRNA in self._mRNAs:
                temp_ribos = copy.copy(mRNA.ribosomes)
                ribosomes_database["mRNA_" + str(mRNA.index).zfill(5)] = [temp_ribos]
            ribosomes_database.close()

    def solve_internal(self, start, end, deltat):
        """
        iterate TRSL for the timerange (start, end, deltat)
        """
        log.info("solve_internal: simulation from %s to %s", start, end)

        fieldnames = ["proteins", "ribos._bound", "ribos._free", "tRNA_bound", "tRNA_free", "ATP", "AMP", "GTP", "GDP", "peptide_bonds"]
        self.initialize_solve_internal(fieldnames)

        self.timerange = np.arange(start, end, deltat)
        for time in self.timerange:
            if int(time) == time:
                print("solve_internal: reached time %.2f sec." % time)
            log.info("################################################################################################")
            log.info(f"solve_internal: time: {time:.2f}")
            log.info("################################################################################################")

            self.update_processes(deltat, time)

            log.info("solve_internal: sum(self.proteins.values()) = %s", sum(self.proteins.values()))
            log.info("solve_internal: protein length:  %s", self.protein_length)
            log.info("solve_internal: bound ribosomes: %s", self.ribo_bound)
            log.info("solve_internal: free ribosomes:  %s", self.ribo_free)
            log.info("solve_internal: collisions/no collisions:  %s/%s", self.collision, self.nocollision)
            log.info("solve_internal: bound tRNA:      %s", sum(self.tRNA_bound.values()))
            log.info("solve_internal: free tRNA:       %s", sum(self.tRNA_free.values()))
            # careful: fieldnames and fieldvalues must be in the same order
            fields = {"proteins": self.proteins, "ribos._bound": self.ribo_bound, "ribos._free": self.ribo_free,
                      "tRNA_bound": sum(self.tRNA_bound.values()), "tRNA_free": sum(self.tRNA_free.values()),
                      "ATP": self.ATP, "AMP": self.AMP, "GTP": self.GTP, "GDP": self.GDP, "peptide_bonds": self.protein_length}
            #fieldvalues = fields.values()

            self.update_solve_internal(deltat, fields, start, time)

    def execute(self, trange, x0 = None):
        """
        :param trange: list of time steps, depending on step size.
                       eg: starttime = 0
                           endtime = 1
                           step size = 0.1
                           --> trange list = [0, 0.1, 0.2, 0.3,....,1.0]
                x0:    initial state vector as a dictionary. It's not given at first call
                       units is number of particles.
                       eg. x0 = {'ATP': 100, 'NADH':5000, 'GLC':700}

        :return  return stateVector, info
                        stateVector: same type as x0 with last particle numbers
                        info:        status message for solverstep
                                     it's a  dictionary like info = {'message': 'everythings fine'}


        :raise NotImplementedError:
        """

        #TODO: implement method
        return None


if __name__ == "__main__":
    log.basicConfig(level=log.INFO, format='%(message)s', stream=sys.stdout)

    trsl = TRSL(nribo=20)
    # overwrite number of transcripts:
    trsl.n_mRNA = 60

    trsl.solve_internal(0.0, 300.0, deltat=0.2)
    trsl.dump_results(description='TRSL_unspecific_low_ribosomes_results')

    '''
    # Profiling:
    cProfile.run('trsl.solve_internal(0.0, 20.0, deltat=1.0)', 'trsl_profile')
    p = pstats.Stats('trsl_profile')
    p.strip_dirs().sort_stats('cumulative').print_stats('TRSL')
    '''


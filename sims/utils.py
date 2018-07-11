import matplotlib.pyplot as plt
import fwdpy11 as fp11
import fwdpy11.wright_fisher as wf
import fwdpy11.model_params as model_params
import numpy as np
import fwdpy11.sampling
import libsequence.polytable as polyt
from libsequence.summstats import PolySIM
import pandas as pd
from libsequence.windows import Windows
import sys
import pickle


def get_demography(path):
    demog = pd.read_csv(path).N
    demog=demog.as_matrix()
    return(demog)


def write_output(recorder,out_path,name,replicate):
    pi = pd.DataFrame(recorder.pi[1:], columns=recorder.pi[0])
    pi['replicate'] = replicate
    pi.to_csv(out_path + '%s_%s_pi.csv' % (replicate,name), index=False)
    singleton = pd.DataFrame(recorder.singleton[1:], columns=recorder.singleton[0])
    singleton['replicate'] = replicate
    singleton.to_csv(out_path + '%s_%s_xi.csv' % (replicate,name), index=False)
    tajimasD= pd.DataFrame(recorder.tajimasD[1:], columns=recorder.tajimasD[0])
    tajimasD['replicate'] = replicate
    tajimasD.to_csv(out_path + '%s_%s_tajD.csv' % (replicate,name), index=False)


def str2byte(tup,fmtstring):
    byte_tup = (tup[0],bytearray(tup[1],fmtstring))
    return(byte_tup)

class neutral_div:
    def __init__(self,set_gen,final,Nstart,nwindows=55,beginning = 50,replicate=1):
        self.nwindows = nwindows
        self.val_per_window = 4
        self.beginning = beginning
        self.pi = [['gen']+[i for i in range((self.nwindows-self.beginning)*self.val_per_window)]]
        self.singleton = [['gen']+[i for i in range((self.nwindows-self.beginning)*self.val_per_window)]]
        self.tajimasD = [['gen']+[i for i in range((self.nwindows-self.beginning)*self.val_per_window)]]
        self.counter = 1
        self.final = final
        self.set_gen = set_gen
        self.Nstart = Nstart
        self.__rng = fp11.GSLrng(np.random.randint(42+float(replicate)))

    def __call__(self, pop):


        if self.counter % 50 == 0  or (pop.generation==self.final):
 #           print(pop.generation)
            actual_gen = np.around([(pop.generation-self.set_gen)/self.Nstart],decimals=3)

            # sample chromosomes from the population
            chr_sampled = 400
            samp = fp11.sampling.sample_separate(self.__rng, pop, chr_sampled, True)
            neutral_sample = polyt.SimData([str2byte(mut, 'utf-8') for mut in samp[0]])

            # split into windows
            w = Windows(neutral_sample, window_size=1/self.val_per_window, step_len=1/self.val_per_window, starting_pos=self.beginning, ending_pos=self.nwindows)

            # calculate summaries
            window_pi = np.around([PolySIM(w[i]).thetapi() for i in range(len(w))],decimals=3)
            window_singleton = np.around([PolySIM(w[i]).numsingletons() for i in range(len(w))],decimals=3)
            window_tajimasD = np.around([PolySIM(w[i]).tajimasd() for i in range(len(w))],decimals=3)

            # add data to output
            self.pi.append(np.append(actual_gen,window_pi))
            self.singleton.append(np.append(actual_gen,window_singleton))
            self.tajimasD.append(np.append(actual_gen,window_tajimasD))
        self.counter += 1


class track_burnin:
    def __init__(self):
        self.counter = 1
    def __call__(self, pop):
        if self.counter % 10000 == 0:
  #          print(pop.generation)
            mut_neut = np.array([(i) for i, j in zip(pop.mcounts, pop.mutations) if
                                 i > 0 and j.neutral is True and j.g == pop.generation and j.pos>0 and j.pos<50])
            mut_sel = np.array([(i) for i, j in zip(pop.mcounts, pop.mutations) if
                               i > 0 and j.neutral is False and j.g == pop.generation and j.pos>0 and j.pos<50])
            s_sel = np.array([(j.s) for i, j in zip(pop.mcounts, pop.mutations) if
                    i > 0 and j.neutral is False and j.g == pop.generation and  j.pos > 0 and j.pos < 50])
            print(pop.generation,mut_sel.sum(),mut_neut.sum(),s_sel.mean())
        self.counter += 1

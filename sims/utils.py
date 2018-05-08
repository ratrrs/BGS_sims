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


def regions_dfe(species='human',neutral=False):
    if species == 'human':
        mu = 1.66e-8
    elif species == 'maize':
        mu = 3.0e-7
    elif species == 'generic':
        mu = 1e-8


    if neutral == False:
        print('with BGS')
        if species == 'human' or species == 'generic':

            sregion = [fp11.GammaS(50, 51, .07, -0.029426,0.184,h=1.0, coupled=True), # coding DFE 7% of all mutations in this region
                      fp11.GammaS(50, 51, 0.13, -0.000518,0.0415,h=1.0, coupled=True) # conserved non-coding DFE 13% of all mutations in this region
                      ]
            nregion = [fp11.Region(i, i + 1, 1, coupled=True) for i in range(50)] + \
                      [fp11.Region(50, 51, 0.8, coupled=True)] + \
                      [fp11.Region(i, i + 1, 1, coupled=True) for i in range(51, 101)]  # 80 % of sites are neutral
            # Mutation rate
            mu_s = mu_n = rec = mu * 20000 * 101
            rates = [mu_s, mu_n, rec]

        elif species =='maize':
            sregion = [fp11.GammaS(10, 11, .2, -0.083, 0.1514, h=1.0, coupled=True)]

            nregion = [fp11.Region(i, i + 1, 1, coupled=True) for i in range(10)] + \
                      [fp11.Region(10, 11, 0.8, coupled=True)] + \
                      [fp11.Region(i, i + 1, 1, coupled=True) for i in range(11, 21)]  # 80 % of sites are neutral

            # Mutation rate
            mu_s = mu_n = rec = mu * 20000 * 21
            rates = [mu_s,mu_n,rec]

    elif neutral== True:
        print('neutral')
        sregion = []
        if species == 'human' or species == 'generic':

            nregion = [fp11.Region(i,i+1,1, coupled=True) for i in range(101)]
            # Mutation rates
            mu_n = rec = mu * 20000 *101
            rates = [mu_n,0,rec]
        elif species == 'maize':
            nregion = [fp11.Region(i,i+1,1, coupled=True) for i in range(21)]
            # Mutation rates
            mu_n = rec = mu * 20000 *21
            rates = [mu_n,0,rec]

    return(sregion,nregion,rates)



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
    def __init__(self,set_gen,final,Nstart,nwindows=101):
        self.nwindows = nwindows
        self.pi = [['gen']+[i for i in range(self.nwindows)]]
        self.singleton = [['gen']+[i for i in range(self.nwindows)]]
        self.tajimasD = [['gen']+[i for i in range(self.nwindows)]]
        self.counter = 1
        self.final = final
        self.set_gen = set_gen
        self.Nstart = Nstart
    def __call__(self, pop):
        rng3 = fp11.GSLrng(np.random.randint(420000))
        if self.counter % 100 == 0  or (pop.generation==self.final):
 #           print(pop.generation)
            actual_gen = np.around([(pop.generation-self.set_gen)/self.Nstart],decimals=3)
            if pop.N < 1000:
            	ind_sampled = pop.N
            else:
            	ind_sampled = 1000
            samp = fp11.sampling.sample_separate(rng3, pop, ind_sampled, True)
            neutral_sample = polyt.SimData([str2byte(mut, 'utf-8') for mut in samp[0]])
            w = Windows(neutral_sample, window_size=1, step_len=1, starting_pos=0., ending_pos=float(self.nwindows))
            window_pi = np.around([PolySIM(w[i]).thetapi() for i in range(len(w))],decimals=3)
            window_singleton = np.around([PolySIM(w[i]).numsingletons() for i in range(len(w))])
            window_tajimasD = np.around([PolySIM(w[i]).tajimasd() for i in range(len(w))],decimals=3)
            self.pi.append(np.append(actual_gen,window_pi))
            self.singleton.append(np.append(actual_gen,window_singleton))
            self.tajimasD.append(np.append(actual_gen,window_tajimasD))
        self.counter += 1


class track_burnin:
    def __init__(self):
        self.counter = 1
    def __call__(self, pop):
        if self.counter % 1000 == 0:
            print(pop.generation)
        self.counter += 1

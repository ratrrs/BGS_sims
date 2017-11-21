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


def regions_human_dfe(mu,neutral=False):
    if neutral == False:
        print('with BGS')
        sregion = [fp11.GammaS(50, 51, .07, -0.029426,0.184,h=0, coupled=True), # coding DFE 7% of all mutations in this region
                  fp11.GammaS(50, 51, 0.13, -0.000518,0.0415,h=0, coupled=True) # conserved non-coding DFE 13% of all mutations in this region
                  ]
        nregion = [fp11.Region(i,i+1,1, coupled=True) for i in range(50)] + \
                  [fp11.Region(50,51,0.8, coupled=True)] +\
                  [fp11.Region(i,i+1,1, coupled=True) for i in range(51,101)]# 80 % of sites are neutral
        # Mutation rates
        mu_s = mu_n = rec = mu * 20000 *101
        rates = [mu_s,mu_n,rec]
    elif neutral== True:
        print('neutral')
        sregion = []
        nregion = [fp11.Region(i,i+1,1, coupled=True) for i in range(101)]
        # Mutation rates
        mu_n = rec = mu * 20000 *101
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
    def __init__(self,set_gen,final,Nstart):
        self.pi = [['gen']+[i for i in range(101)]]
        self.singleton = [['gen']+[i for i in range(101)]]
        self.tajimasD = [['gen']+[i for i in range(101)]]
        self.counter = 1
        self.final = final
        self.set_gen = set_gen
        self.Nstart = Nstart
    def __call__(self, pop):
        rng3 = fp11.GSLrng(np.random.randint(420000))
        if self.counter % 100 == 0  or (pop.generation==self.final):
            print(pop.generation)
            samp = fp11.sampling.sample_separate(rng3, pop, 1000, True)
            neutral_sample = polyt.SimData([str2byte(mut, 'utf-8') for mut in samp[0]])
            w = Windows(neutral_sample, window_size=1, step_len=1, starting_pos=0., ending_pos=101.0)
            window_pi = np.around([PolySIM(w[i]).thetapi() for i in range(len(w))],decimals=3)
            window_singleton = np.around([PolySIM(w[i]).numsingletons() for i in range(len(w))])
            window_tajimasD = np.around([PolySIM(w[i]).tajimasd() for i in range(len(w))],decimals=3)

            self.pi.append(np.append([(pop.generation-self.set_gen)/self.Nstart],window_pi))
            self.singleton.append(np.append([(pop.generation-self.set_gen)/self.Nstart],window_singleton))
            self.tajimasD.append(np.append([pop.generation-self.set_gen/self.Nstart],window_tajimasD))
        self.counter += 1



def regions_human_dfe_torres(mu,neutral=False):
    if neutral == False:
        print('with BGS')
        sregion = [fp11.GammaS(i, i+1, .07, -0.029426,0.184, coupled=True) for i in range(50)]+ \
                  [fp11.GammaS(i, i+1, 0.13, -0.000518, 0.0415, coupled=True) for i in range(50)] + \
                  [fp11.GammaS(i, i+1, .07, -0.029426, 0.184, coupled=True) for i in range(51,101)] + \
                  [fp11.GammaS(i, i+1, 0.13, -0.000518, 0.0415, coupled=True) for i in range(51,101)]

        nregion = [fp11.Region(i,i+1,0.8, coupled=True) for i in range(50)] + \
                  [fp11.Region(50,51,1, coupled=True)] +\
                  [fp11.Region(i,i+1,0.8, coupled=True) for i in range(51,101)]# 80 % of sites are neutral
        # Mutation rates
        mu_s = mu_n = rec = mu * 20000 *101
        rates = [mu_s,mu_n,rec]
    elif neutral== True:
        print('neutral')
        sregion = []
        nregion = [fp11.Region(i,i+1,1, coupled=True) for i in range(101)]
        # Mutation rates
        mu_n = rec = mu * 20000 *101
        rates = [mu_n,0,rec]
    return(sregion,nregion,rates)
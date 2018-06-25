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
    elif species == 'generic' or 'test':
        mu = 3.0e-8


    if neutral == False:
        print('with BGS')
        if species == 'human':
            sregion = [fp11.GammaS(10, 11, 1.,-0.029426,0.184,h=1.0, coupled=True), # coding DFE 1/3 of sel mutations in this region
                      fp11.GammaS(10, 11, 2.,-0.000518,0.0415,h=1.0, coupled=True) # conserved non-coding DFE 2/3 of sel mutations in this region
                      ]
            nregion = [fp11.Region(i, i + 1, 1., coupled=True) for i in range(10)] + \
                      [fp11.Region(10, 11, 0.20, coupled=True)] + \
                      [fp11.Region(i, i + 1, 1., coupled=True) for i in range(11, 21)]  # 80 % of sites are neutral
            print('Number of neutral regions:',len(nregion))

            # Mutation rate
            rec = mu * 20000 * 21
            mu_s = mu * 20000 *.8
            mu_n = mu * 20000 * 20 + mu * 20000 * 0.8
            rates = [mu_n, mu_s, rec]

        elif species =='maize':
            sregion = [fp11.GammaS(10, 11, 1, -0.83, 0.1514, h=1.0, coupled=True)]
            nregion = [fp11.Region(i, i + 1, 1, coupled=True) for i in range(10)] + \
                      [fp11.Region(10, 11, .2, coupled=True)] + \
                      [fp11.Region(i, i + 1, 1, coupled=True) for i in range(11, 21)]  # 20 % of sites are neutral
            print('Number of neutral regions:',len(nregion))
            # Mutation rate
            mu_s = mu * 20000 *0.8
            mu_n = mu * 20000 *20 + mu * 20000 * 0.2
            rec = mu * 20000 * 21
            rates = [mu_n,mu_s,rec]

        elif species== 'test' or species == 'generic':
            sregion = [fp11.GammaS(10, 11, 1, -0.83, 0.01514, h=1.0, coupled=True)]
            nregion = [fp11.Region(i, i + 1, 1., coupled=True) for i in range(10)] + \
                      [fp11.Region(10, 11, 0.2, coupled=True)] + \
                      [fp11.Region(i, i + 1, 1., coupled=True) for i in range(11, 21)]  # 20 % of sites are neutral
            print('Number of neutral regions:', len(nregion))
            # Mutation rate
            rec = mu * 20000 * 21
            mu_s = mu * 20000 *.8
            mu_n = mu * 20000 * 20 + mu * 20000 * 0.2
            rates = [mu_n, mu_s, rec]

    elif neutral== True:
        print('neutral')
        sregion = []
        if species == 'human':
            nregion = [fp11.Region(i,i+1,1, coupled=True) for i in range(21)]
            # Mutation rates
            mu_n = rec = mu * 20000 *21
            rates = [mu_n,0,rec]

        elif species == 'maize' or species == 'test' or species == 'generic':
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
    def __init__(self,set_gen,final,Nstart,nwindows=21):
        self.nwindows = nwindows
        self.val_per_window = 10
        self.pi = [['gen']+[i for i in range(self.nwindows*self.val_per_window)]]
        self.singleton = [['gen']+[i for i in range(self.nwindows*self.val_per_window)]]
        self.tajimasD = [['gen']+[i for i in range(self.nwindows*self.val_per_window)]]
        self.counter = 1
        self.final = final
        self.set_gen = set_gen
        self.Nstart = Nstart
    def __call__(self, pop):
        rng3 = fp11.GSLrng(np.random.randint(420000))
        if self.counter % 50 == 0  or (pop.generation==self.final):
 #           print(pop.generation)
            actual_gen = np.around([(pop.generation-self.set_gen)/self.Nstart],decimals=3)
            ind_sampled = 400
            samp = fp11.sampling.sample_separate(rng3, pop, ind_sampled, True)
            neutral_sample = polyt.SimData([str2byte(mut, 'utf-8') for mut in samp[0]])
            w = Windows(neutral_sample, window_size=1/self.val_per_window, step_len=1/self.val_per_window, starting_pos=0., ending_pos=float(self.nwindows))
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
        if self.counter % 10000 == 0:
  #          print(pop.generation)
            mut_neut = np.array([(i) for i, j in zip(pop.mcounts, pop.mutations) if
                                 i > 0 and j.neutral is True and j.g == pop.generation and j.pos>10 and j.pos<11])
            mut_sel = np.array([(i) for i, j in zip(pop.mcounts, pop.mutations) if
                               i > 0 and j.neutral is False and j.g == pop.generation and j.pos>10 and j.pos<11])
            s_sel = np.array([(j.s) for i, j in zip(pop.mcounts, pop.mutations) if
                    i > 0 and j.neutral is False and j.g == pop.generation and  j.pos > 10 and j.pos < 11])
            print(pop.generation,mut_sel.sum(),mut_neut.sum(),s_sel.mean())
        self.counter += 1

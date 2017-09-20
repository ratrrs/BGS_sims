#!/bin/env python

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

def str2byte(tup,fmtstring):
    byte_tup = (tup[0],bytearray(tup[1],fmtstring))
    return(byte_tup)

class neutral_div:
    def __init__(self):
        self.values = []#[['gen']+[i for i in range(25)]]
        self.singleton = []#[['gen'] + [i for i in range(25)]]
    def __call__(self, pop):
#        results = np.array([i.s for i in pop.mutations])
#        print(results)
        if pop.generation % 100 == 0:
            samp = fp11.sampling.sample_separate(rng2, pop, 500, True)
            neutral_sample = polyt.SimData([str2byte(mut, 'utf-8') for mut in samp[0]])
            w = Windows(neutral_sample, window_size=1, step_len=1, starting_pos=0., ending_pos=25.0)
            window_pi = [PolySIM(w[i]).thetapi() for i in range(len(w))]
            window_singleton = [PolySIM(w[i]).numsingletons() for i in range(len(w))]
            self.values.append([pop.generation]+window_pi)
            self.singleton.append([pop.generation]+window_singleton)
#start = int(sys.argv[1])
#end = int(sys.argv[2])
#print(start,end)
#for i in range(start,end):

run = int(sys.argv[1])

# European Human Demography
Nstart = int(1e4)
OOA_bottleneck_size = int(0.322*Nstart)
OOA_time = int(0.109*2*Nstart) # 2180
European_bottleneck_size = 0.185*Nstart
European_bottleneck_time = int(0.053*2*Nstart) #1060
final_European_size = 4.44*Nstart

before = np.array([Nstart]*200,dtype=np.uint32)
OOA_bottleneck = np.array([OOA_bottleneck_size]*OOA_time,dtype=np.uint32)
x = np.linspace(np.log(European_bottleneck_size), np.log(final_European_size), European_bottleneck_time)
European_bottleneck = np.exp(x).round().astype(np.int32)
demog_euro = np.concatenate((before,OOA_bottleneck,European_bottleneck)).astype(np.int32)
print(len(demog_euro))
# Maize demography

Npresent = 3*Nstart
Nbneck = 0.05*Nstart
x = np.linspace(np.log(Nbneck), np.log(Npresent), 1000)
bneckpop = np.exp(x).round().astype(np.int32)

before = np.array([Nstart]*2440,dtype=np.uint32) # Matched for Human demography
demog_maize = np.concatenate((before,bneckpop)).astype(np.int32)
print(len(demog_maize))
# burnin demography
nlist=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)


rng2=fp11.GSLrng(42*run)
pop=fp11.SlocusPop(Nstart)


# regions
#sregion = [fp11.GammaS(12, 13, 0.75, -5e-4,0.3, coupled=True)]
sregion = []
nregion = [fp11.Region(i,i+1,1, coupled=True) for i in range(25)]

p = {'nregions':nregion,
'sregions': sregion,
'recregions':[fp11.Region(0,25,1)],
'rates':(0.03,0,0.03),
'demography':nlist,
}

params = fp11.model_params.SlocusParams(**p)

# run burnin

wf.evolve(rng2, pop,params)

ppop = pickle.dumps(pop,-1)
#Unpickle to create a new pop:
pop2 = pickle.loads(ppop)
print(pop==pop2)

print('burnin done')

# run Maize
p['demography'] = demog_maize
params = fp11.model_params.SlocusParams(**p)
rec1=neutral_div()
print('GEN maize',pop2.generation)

wf.evolve(rng2, pop,params,rec1)

output = pd.DataFrame(rec1.values[1:],columns=rec1.values[0])
output['run'] = run
output.to_csv("results/tmp_neutral2/sim_neutralPi_maize_all_neutral_%d.csv"%run,index=False,header=False)


singleton = pd.DataFrame(rec1.singleton[1:],columns=rec1.singleton[0])
singleton['run'] = run
singleton.to_csv("results/tmp_neutral2/sim_singleton_maize_all_neutral_%d.csv"%run,index=False,header=False)


# run Europenan humans
p['demography'] = demog_euro
params = fp11.model_params.SlocusParams(**p)
rec1=neutral_div()
print('GEN euro',pop2.generation)
wf.evolve(rng2, pop2,params,rec1)

output = pd.DataFrame(rec1.values[1:],columns=rec1.values[0])
output['run'] = run
output.to_csv("results/tmp_neutral2/sim_neutralPi_European_all_neutral_%d.csv"%run,index=False,header=False)


singleton = pd.DataFrame(rec1.singleton[1:],columns=rec1.singleton[0])
singleton['run'] = run
singleton.to_csv("results/tmp_neutral2/sim_singleton_European_all_neutral_%d.csv"%run,index=False,header=False)


print('finished run %d'%run)

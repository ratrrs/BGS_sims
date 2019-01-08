from utils import *
import pandas as pd
import sys
import os
import lzma
import fwdpy11 as fp11

Nstart = int(sys.argv[1]) # 400 
out_path = sys.argv[2]#'../results/fixedn/n400'
replicate = str(sys.argv[3]) #replicate name

print(Nstart)
mu = 1.66e-8

class div_rec:
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


        if pop.generation==self.final:
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


################## simulate neutral ############################
recregion =[fp11.Region(50,55,1., coupled=True)]
sregion= []
nregion = [fp11.Region(50, 55, 1., coupled=True)]

# Mutation rates
mu_n = mu * 40000 * 5
rec = 8.2e-10 * 40000 * 5
rates = [mu_n, 0, rec]

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*float(Nstart)),dtype=np.uint32)

mypop =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(122*(1+float(replicate))))

print('rec')
[print(i) for i in recregion]
print('sregion')
[print(i) for i in sregion]
print('nregion')
[print(i) for i in nregion]
print(rates)

p = {'nregions':nregion,
'sregions': sregion,
'recregions':recregion,
'rates':rates,
'demography':burnin,
}
params = fp11.model_params.SlocusParams(**p)

rec1 = div_rec(0, final=int(10*Nstart), Nstart=Nstart,replicate=replicate)



# simulate until equilibrium
wf.evolve(rng2, mypop,params,rec1)
write_output(rec1, model_path, 'neutral', replicate)




#ppop = pickle.dumps(mypop,-1)
# pickle equilibirum population
#burnin_name = out_path + "burnins/burnin_neut_%s.lzma9" % replicate
#with lzma.open(burnin_name, "wb", preset=9) as f:
#    pickle.dump(mypop, f, -1)

print('burnin done')
print('Generation',mypop.generation)


################## simulate BGS ############################

recregion =[fp11.Region(0,55,1., coupled=True)]
sregion =   [fp11.GammaS(0, 50, 1., -0.029426, 0.184, h=1.0, coupled=True)] +\
            [fp11.GammaS(0, 50, 2., -0.000518, 0.0415, h=1.0, coupled=True)] # conserved non-coding DFE 2/3 of sel mutations in this region

#sregion = [fp11.GammaS(-1, 0, 1, -0.83, 0.01514, h=1.0, coupled=True)]


nregion = [fp11.Region(50, 55, 1., coupled=True)]

# Mutation rate
rec = 8.2e-10 * 40000 * 55
mu_s = mu * 40000 * 50 * 0.2
mu_n = mu * 40000 * 5
rates = [mu_n, mu_s, rec]

# constant size for 10 N generations
#burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

mypop =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(42*(1+float(replicate))))


p = {'nregions':nregion,
'sregions': sregion,
'recregions':recregion,
'rates':rates,
'demography':burnin,
}

params = fp11.model_params.SlocusParams(**p)

burn_rec = track_burnin()
print('sregions')
[print(i) for i in p['sregions']]
print('nregions')
[print(i) for i in p['nregions']]
print('Recomb region')
[print(i) for i in p['recregions']]
print(p['rates'])

# simulate until equilibrium

rec1 = div_rec(0, final=int(10*Nstart), Nstart=Nstart,replicate=replicate)


wf.evolve(rng2, mypop,params,rec1)
print('burnin done')
print('Generation',mypop.generation)

write_output(rec1, model_path, 'bgs', replicate)
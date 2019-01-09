from utils import *
import pandas as pd
import sys
import os
import lzma
import fwdpy11 as fp11

out_path = sys.argv[1]#'../results/fixedn/'
replicate = str(sys.argv[2]) #replicate name

starts = [200,400,800,1000,1500,2500,5000,10000,15000,20000,25000,30000]

mu = 1.66e-8


def write_output(recorder,out_path,name,replicate):
    pi = pd.DataFrame(recorder.pi[1:], columns=recorder.pi[0])
    pi['replicate'] = replicate
    pi_file =out_path + '%s_%s_pi.csv' % (replicate,name)
    if not os.path.isfile(pi_file):
    	pi.to_csv(pi_file, index=False)
    else:
    	pi.to_csv(pi_file,mode='a',header=False, index=False)
    singleton = pd.DataFrame(recorder.singleton[1:], columns=recorder.singleton[0])
    singleton['replicate'] = replicate
    xi_file = out_path + '%s_%s_xi.csv' % (replicate,name)
    if not os.path.isfile(xi_file):
    	singleton.to_csv(xi_file, index=False)
    else:
    	singleton.to_csv(xi_file,mode='a',header=False, index=False)
#    tajimasD= pd.DataFrame(recorder.tajimasD[1:], columns=recorder.tajimasD[0])
#    tajimasD['replicate'] = replicate
#    tajimasD.to_csv(out_path + '%s_%s_tajD.csv' % (replicate,name), index=False)




class div_rec:
    def __init__(self,set_gen,final,Nstart,nwindows=55,beginning = 50,replicate=1):
        self.nwindows = nwindows
        self.val_per_window = 4
        self.beginning = beginning
        self.pi = [['N']+[i for i in range((self.nwindows-self.beginning)*self.val_per_window)]]
        self.singleton = [['N']+[i for i in range((self.nwindows-self.beginning)*self.val_per_window)]]
        self.tajimasD = [['N']+[i for i in range((self.nwindows-self.beginning)*self.val_per_window)]]
        self.counter = 1
        self.final = final
        self.set_gen = set_gen
        self.Nstart = Nstart
        self.__rng = fp11.GSLrng(np.random.randint(42+float(replicate)))

    def __call__(self, pop):


        if pop.generation==10*Nstart:
            print("record at:",pop.generation)
            actual_gen = Nstart

            # sample chromosomes from the population
            chr_sampled = 400
            samp = fp11.sampling.sample_separate(self.__rng, pop, chr_sampled, True)
            neutral_sample = polyt.SimData([str2byte(mut, 'utf-8') for mut in samp[0]])

            # split into windows
            w = Windows(neutral_sample, window_size=1/self.val_per_window, step_len=1/self.val_per_window, starting_pos=self.beginning, ending_pos=self.nwindows)

            # calculate summaries
            window_pi = np.around([PolySIM(w[i]).thetapi() for i in range(len(w))],decimals=3)
            window_singleton = np.around([PolySIM(w[i]).numsingletons() for i in range(len(w))],decimals=3)
#            window_tajimasD = np.around([PolySIM(w[i]).tajimasd() for i in range(len(w))],decimals=3)

            # add data to output
            self.pi.append(np.append(actual_gen,window_pi))
            self.singleton.append(np.append(actual_gen,window_singleton))
#            self.tajimasD.append(np.append(actual_gen,window_tajimasD))
        self.counter += 1


################## simulate neutral ############################
recregion =[fp11.Region(50,55,1., coupled=True)]
sregion= []
nregion = [fp11.Region(50, 55, 1., coupled=True)]

# Mutation rates
mu_n = mu * 40000 * 5
rec = 8.2e-10 * 40000 * 5
rates = [mu_n, 0, rec]

for i in starts:
	Nstart = i
	print(Nstart)
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
	write_output(rec1, out_path, 'neutral', replicate)




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

for i in starts:
	Nstart=i
	print(Nstart)
	mypop =  fp11.SlocusPop(Nstart)

	#prepare random number gernerator
	rng2 = fp11.GSLrng(np.random.randint(122*(1+float(replicate))))


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

	write_output(rec1, out_path, 'bgs', replicate)
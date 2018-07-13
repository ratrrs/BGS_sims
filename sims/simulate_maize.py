from utils import *
import pandas as pd
import sys
import os
import lzma
import fwdpy11 as fp11

demog_file = sys.argv[1]#'../demographies/generic_models.csv'
out_path = sys.argv[2]#'../results/generic/'
replicate = str(sys.argv[3]) #replicate name


#demographies = pd.read_csv(demog_file)
#models = np.array(demographies.columns[:-1])
#print(models)

demog =get_demography(demog_file)


Nstart = int(demog.item(0))
mu = 1.66e-7




################## simulate neutral ############################
recregion =[fp11.Region(50,55,1., coupled=True)]
sregion= []
nregion = [fp11.Region(50, 55, 1., coupled=True)]
# Mutation rates
mu_n = mu * 40000 * 5
rec = 8.2e-9 * 40000 * 5
rates = [mu_n, 0, rec]

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

mypop =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(42*(1+float(replicate))))

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

burn_rec = track_burnin()
# simulate until equilibrium
wf.evolve(rng2, mypop,params,burn_rec)




#ppop = pickle.dumps(mypop,-1)
# pickle equilibirum population
burnin_name = out_path + "burnins/burnin_neut_%s.lzma9" % replicate
with lzma.open(burnin_name, "wb", preset=9) as f:
    pickle.dump(mypop, f, -1)

print('burnin done')
print('Generation',mypop.generation)




with lzma.open(burnin_name, 'rb') as f:
    pop2 = pickle.load(f)
print(mypop==pop2)
print(pop2.generation)

p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10 * Nstart) + 200  # adjust generation labels without burnin and start
rec1 = neutral_div(set_gen, final=pop2.generation + len(demog) + 200, Nstart=Nstart)

wf.evolve(rng2, pop2, params, rec1)
print('Generation', pop2.generation)

# write output
write_output(rec1, out_path, 'neutral', replicate)



################## simulate BGS ############################

recregion =[fp11.Region(0,55,1., coupled=True)]
sregion =   [fp11.GammaS(0, 50, 1., -0.29426, 0.184, h=1.0, coupled=True)] +\
            [fp11.GammaS(0, 50, 2., -0.00518, 0.0415, h=1.0, coupled=True)] # conserved non-coding DFE 2/3 of sel mutations in this region

#sregion = [fp11.GammaS(-1, 0, 1, -0.83, 0.01514, h=1.0, coupled=True)]


nregion = [fp11.Region(50, 55, 1., coupled=True)]

# Mutation rate
rec = 8.2e-9 * 40000 * 55
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


wf.evolve(rng2, mypop,params,burn_rec)
print('burnin done')
print('Generation',mypop.generation)


#ppop = pickle.dumps(mypop,-1)
# pickle equilibirum population
burnin_name = out_path + "burnins/burnin_bgs_%s.lzma9" % replicate
with lzma.open(burnin_name, "wb", preset=9) as f:
    pickle.dump(mypop, f, -1)




with lzma.open(burnin_name, 'rb') as f:
    pop2 = pickle.load(f)
print(pop2.generation)
#Unpickle to create a new pop:
#    pop2 = pickle.loads(ppop)
#    print(mypop==pop2)

p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10 * Nstart) + 200  # adjust generation labels without burnin and start
rec1 = neutral_div(set_gen, final=pop2.generation + len(demog) + 200, Nstart=Nstart)

wf.evolve(rng2, pop2, params, rec1)
print('Generation', pop2.generation)

# write output
write_output(rec1, out_path, 'bgs', replicate)
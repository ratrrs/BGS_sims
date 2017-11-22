from utils import *
import pandas as pd
import sys

demog_file = sys.argv[1]#'../demographies/tennessen.csv'
out_path = sys.argv[2]#'results/'
replicate = str(sys.argv[3])
# get the demography

demog =get_demography(demog_file)


################## simulate neutral ############################
sregion,nregion,rates = regions_human_dfe(1.66e-8,neutral = True)
# define starting population size (get first
Nstart = int(demog.item(0))

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

humans =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(420000))


p = {'nregions':nregion,
'sregions': sregion,
'recregions':[fp11.Region(0,len(nregion),1)],
'rates':rates,
'demography':burnin,
}
params = fp11.model_params.SlocusParams(**p)

# simulate until equilibrium
wf.evolve(rng2, humans,params)
print('burnin done')
print('Generation',humans.generation)


# simulate with demography from 10 N generations
p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10*Nstart)+200 # adjust generation labels without burnin and start
rec1=neutral_div(set_gen,final=humans.generation+len(demog)+200,Nstart=Nstart)


wf.evolve(rng2, humans,params,rec1)
print('Generation',humans.generation)

# write output
write_output(rec1,out_path,'neutral',replicate)

################## simulate BGS ############################

sregion,nregion,rates = regions_human_dfe(1.66e-8,neutral = False)



# define starting population size (get first
Nstart = int(demog.item(0))

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

humans =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(420000))


p = {'nregions':nregion,
'sregions': sregion,
'recregions':[fp11.Region(0,len(nregion),1)],
'rates':rates,
'demography':burnin,
}
params = fp11.model_params.SlocusParams(**p)
print('Generation start',humans.generation)

# simulate until equilibrium
wf.evolve(rng2, humans,params)
print('burnin done')
print('Generation',humans.generation)


# simulate with demography from 10 N generations
p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10*Nstart)+200 # adjust generation labels without burnin and start
rec1=neutral_div(set_gen,final=humans.generation+len(demog)+200,Nstart=Nstart)


wf.evolve(rng2, humans,params,rec1)
print('Generation',humans.generation)

# write output
write_output(rec1,out_path,'bgs',replicate)

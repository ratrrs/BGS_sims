from utils import *
import pandas as pd
import sys

demog_file = sys.argv[1]#'../demographies/tennessen.csv'
out_path = sys.argv[2]#'results/'
species = str(sys.argv[3]) #human, maize or generic
replicate = str(sys.argv[4]) #replicate name

# get the demography

demog =get_demography(demog_file)

# define starting population size (get first value from demog file)
Nstart = int(demog.item(0))


print('Simulating neutral and BGS model with following parameters')
print('Species:',species)
print('Ancestral population size:',int(demog.item(0)))
print('total generations', len(demog))


################## simulate neutral ############################

sregion,nregion,rates = regions_dfe(species=species,neutral = True)

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

mypop =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(420000))


p = {'nregions':nregion,
'sregions': sregion,
'recregions':[fp11.Region(0,len(nregion),1)],
'rates':rates,
'demography':burnin,
}
params = fp11.model_params.SlocusParams(**p)

burn_rec = track_burnin()


# simulate until equilibrium
wf.evolve(rng2, mypop,params,burn_rec)
print('burnin done')
print('Generation',mypop.generation)


# simulate with demography from 10 N generations
p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10*Nstart)+200 # adjust generation labels without burnin and start
rec1=neutral_div(set_gen,final=mypop.generation+len(demog)+200,Nstart=Nstart)


wf.evolve(rng2, mypop,params,rec1)
print('Generation',mypop.generation)

# write output
write_output(rec1,out_path,'neutral',replicate)

################## simulate BGS ############################

sregion,nregion,rates = regions_dfe(species=species,neutral = False)



# define starting population size (get first
Nstart = int(demog.item(0))

# constant size for 10 N generations
burnin=np.array([Nstart]*int(10*Nstart),dtype=np.uint32)

mypop =  fp11.SlocusPop(Nstart)

#prepare random number gernerator
rng2 = fp11.GSLrng(np.random.randint(420000))


p = {'nregions':nregion,
'sregions': sregion,
'recregions':[fp11.Region(0,len(nregion),1)],
'rates':rates,
'demography':burnin,
}
params = fp11.model_params.SlocusParams(**p)
print('Generation start',mypop.generation)

# simulate until equilibrium
wf.evolve(rng2, mypop,params)
print('burnin done')
print('Generation',mypop.generation)


# simulate with demography from 10 N generations
p['demography'] = demog
params = fp11.model_params.SlocusParams(**p)

# add recorder that records pi, singletons and tajimas D
set_gen = (10*Nstart)+200 # adjust generation labels without burnin and start
rec1=neutral_div(set_gen,final=mypop.generation+len(demog)+200,Nstart=Nstart)


wf.evolve(rng2, mypop,params,rec1)
print('Generation',mypop.generation)

# write output
write_output(rec1,out_path,'bgs',replicate)

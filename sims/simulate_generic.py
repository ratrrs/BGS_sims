from utils import *
import pandas as pd
import sys
import os
import lzma

demog_file = sys.argv[1]#'../demographies/generic_models.csv'
out_path = sys.argv[2]#'../results/generic/'
replicate = str(sys.argv[3]) #replicate name


demographies = pd.read_csv(demog_file)
models = np.array(demographies.columns[:-1])
#print(models)
Nstart = 20000


################## simulate neutral ############################

sregion,nregion,rates = regions_dfe(species='generic',neutral=True)

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




#ppop = pickle.dumps(mypop,-1)
# pickle equilibirum population
burnin_name = out_path + "burnins/burnin_neut_%s.lzma9" % replicate
with lzma.open(burnin_name, "wb", preset=9) as f:
    pickle.dump(mypop, f, -1)

print('burnin done')
print('Generation',mypop.generation)

for model in models:
    print(model)
    model_id = model.replace(" ", "").lower() +'/'
    print(model_id)
    model_path = out_path + model_id
    if not os.path.exists(model_path):
        os.makedirs(model_path)
    demog = demographies[model].as_matrix()
    demog = demog[~np.isnan(demog)]

    #Unpickle to create a new pop:
    #pop2 = pickle.loads(ppop)

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
    write_output(rec1, model_path, 'neutral', replicate)



################## simulate BGS ############################

sregion,nregion,rates = regions_dfe(species='generic',neutral=False)

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


#ppop = pickle.dumps(mypop,-1)
# pickle equilibirum population
burnin_name = out_path + "burnins/burnin_bgs_%s.lzma9" % replicate
with lzma.open(burnin_name, "wb", preset=9) as f:
    pickle.dump(mypop, f, -1)


for model in models:
    print('simulating BGS in %s'%model)
    model_id = model.replace(" ", "").lower() +'/'
    model_path = out_path + model_id
    demog = demographies[model].as_matrix()
    demog = demog[~np.isnan(demog)]

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
    write_output(rec1, model_path, 'bgs', replicate)
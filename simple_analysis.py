import msprime
import pyslim
import os
import tskit
import allel
import numpy as np
import pandas as pd
import itertools
import glike

path="/Users/olj5016/Documents/arg_selection/"

os.chdir(path)

scaling=10 # scale parameters by 10
Ne=10000 # unscaled Ne
rr=1e-8 # unscaled recombination rate
mr=1e-8 # unscaled mutation rate
s=0.06 # selection coefficient (0.0 for neutral)
sampleSize=8 # number of indidivuals remebered in tree per generation
selPop=2 # subset value for subpopulation in SLiM 2 for p2, 4 for p22
selTime=17500 # time of selection (generations)
selEnd=20000 # time of selection (generations)
cF = 0.1 # condntional frequency of selected allele (only active when s>0.0)
admixture=0 ## admixture proportion set to 0 to turn admixture off

# sN=Ne/scaling # scaled N
# sR=rr*scaling # scaled rr
# sM=mr*scaling # scaled mr
# sS=s*scaling # scaled s

params="s{0}_sT{1}_sE{2}_sP{3}_cF{4}_admix{5}".format(s,selTime, selEnd,selPop, cF, admixture)

burnin = msprime.sim_ancestry(samples=Ne, population_size=Ne, recombination_rate=rr, sequence_length=1e7)
burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
burnin_ts.dump("{0}burnin_simple_{1}.trees".format(path,params))

cmd = "slim -d s=" + str(s) + " -d admix="+ str(admixture)+" -d sampleSize=" + str(sampleSize)+ " -d selPop=" + str(selPop)+ " -d selTime=" + str(selTime) + " -d cF=" + str(cF)+ " ~/arg_selection/simple_gross.slim"
print(cmd)
os.system(cmd)


ts=tskit.load("simplegross_{0}.trees".format(params))

## SUMMARISE INDIVIDUALS - obtain metadata for individuals 'remembered' in ts
rows_list = []
        # run through individuals (inds) in ts
for ind in ts.individuals():
    #print(ind)
    dict1 = {}
            # individual's id in ts
    dict1.update({"id" : ind.id})

    dict1.update({"pop" : ind.population })

    dict1.update({"time" : ind.time})

        #     # genotypes individuals contained (2 nodes, each directing to a haploid genotype)
    dict1.update({"nodes": ind.nodes})

    rows_list.append(dict1)
        # convert from dictionary to data frame (fast)
ind_met = pd.DataFrame(rows_list)

ind_times = np.unique(ts.individual_times).astype(int)







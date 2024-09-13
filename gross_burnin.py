import msprime
import pyslim
import os
import tskit
import allel
import numpy as np
import pandas as pd
import itertools

scaling=10 # scale parameters by 10
Ne=10000 # unscaled Ne
rr=1e-8 # unscaled recombination rate
mr=1e-8 # unscaled mutation rate
s=0.1 # selection coefficient (0.0 for neutral)
sampleSize=4 # number of indidivuals remebered in tree per generation
selPop=2 # subset value for subpopulation in SLiM 2 for p2, 4 for p22
selTime=500 # time of selection (generations)
cF = 0.1 # condntional frequency of selected allele (only active when s>0.0)

sN=Ne/scaling # scaled N
sR=rr*scaling # scaled rr
sM=mr*scaling # scaled mr
sS=s*scaling # scaled s

burnin = msprime.sim_ancestry(samples=sN, population_size=sN, recombination_rate=sR, sequence_length=1e7)
burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
burnin_ts.dump("/Users/olj5016/Documents/arg_selection/gross_burnin.trees")

params="p{0}_t{1}_s{2}_f{3}_sS{4}.txt".format(selPop, selTime, sS, cF, sampleSize)

cmd = "slim -d s=" + str(sS) + " -d sampleSize=" + str(sampleSize)+ " -d selPop=" + str(selPop)+ " -d selTime=" + str(selTime) + " -d cF=" + str(cF)+ " ~/arg_selection/gross_demography.slim"
print(cmd)
os.system(cmd)

ts = tskit.load("/Users/olj5016/Documents/arg_selection/simtree_{0}.trees".format(params))

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

mut_ts = msprime.sim_mutations(ts, rate=sM, discrete_genome=True, keep=True)
rows_list2 = []
for p in set(ind_met['pop']):
     
     for t in np.unique(ind_met[ind_met["pop"]==p]['time']):
            
            sample_ind = pyslim.individuals_alive_at(mut_ts, t, population=p)
            sample = ind_met["nodes"].iloc[sample_ind]

            ## convert to list
            samples= list(itertools.chain(*sample))

            ## create ts of just samples at time t using list
            samp_ts = mut_ts.simplify(samples = samples)

            ## create genotype array to put into scikit-allel  shape ==(var, ind)
            samp_gm=samp_ts.genotype_matrix()
            h= allel.HaplotypeArray(samp_gm)
            samp_ac = h.count_alleles()
            sites=samp_ts.sites_position
            for i in range(len(sites)):
                dict2 = {}
                dict2.update({"pop": p, "time":pyslim.slim_time(mut_ts, t), "pos":sites[i], "A":samp_ac[i,0], "a":samp_ac[i,1]})
                rows_list2.append(dict2)
alleleCounts=pd.DataFrame(rows_list2)
filename="/Users/olj5016/Documents/arg_selection/ac_{0}.txt".format(params)
alleleCounts.to_string(buf = filename, index=False)

reformat="Rscript /Users/olj5016/arg_selection/gross_analysis.R {0}".format(filename)
os.system(reformat)

grosscmd="Rscript /Users/olj5016/Documents/arg_selection/GRoSS-master/GRoSS.R -e {0} -d /Users/olj5016/Documents/arg_selection/treeGross.dot -o /Users/olj5016/Documents/arg_selection/gross_out_{1}.tsv".format(filename, params)
os.system(grosscmd)

    
#from IPython.display import display
#x_limits = [4950000, 5050000]
#ts.draw_svg(path="/Users/olj5016/Documents/arg_selection/gross_arg.pdf",y_axis=True, x_lim=x_limits)



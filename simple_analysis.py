import msprime
import pyslim
import os
import sys
import tskit
#import allel
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import statistics

path="/Users/olj5016/Documents/arg_selection/"
#path="/Users/olivia/Documents/arg_selection/"


os.chdir(path)

sys.path.insert(1, "/Users/olj5016/arg_selection/")
#sys.path.insert(1, "/Users/olivia/arg_selection/")
import glike
import estimate
import miscellaneous


Ne=10000 # unscaled Ne
rr=1e-8 # unscaled recombination rate
mr=1e-8 # unscaled mutation rate
s=1.0 # selection coefficient (0.0 for neutral)
sampleSize=8 # number of indidivuals remebered in tree per generation
selPop=2 # subset value for subpopulation in SLiM 2 for p2, 4 for p22
selTime=17500 # time of selection (generations)
selEnd=20000 # time of selection (generations)
cF = 0.5 # conditional frequency of selected allele (only active when s>0.0)
admixture=0.0 ## admixture proportion set to 0 to turn admixture off
rep=0


params="{6}_s{0}_sT{1}_sE{2}_sP{3}_cF{4}_admix{5}".format(s,selTime, selEnd,selPop, cF, admixture,rep)

if os.path.isfile("{0}simplegross_{1}.trees".format(path,params))==False:
    
    burnin = msprime.sim_ancestry(samples=Ne, population_size=Ne, recombination_rate=rr, sequence_length=1e7)
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1, stage="late")
    burnin_ts.dump("{0}burnin_simple_{1}.trees".format(path,params))
    
    cmd = "slim -d s=" + str(s) + " -d rep="+str(rep)+" -d admix="+ str(admixture)+" -d sampleSize=" + str(sampleSize)+ " -d selPop=" + str(selPop)+ " -d selTime=" + str(selTime) + " -d cF=" + str(cF)+ " ~/arg_selection/simple_gross.slim"
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

mut_ts = msprime.sim_mutations(ts, rate=mr, discrete_genome=True, keep=True)

modernInds=ts.samples(time=1) ## sample inds from final generation
mp2=ts.samples(2,time=1)

mod_ts=mut_ts.simplify(samples=mp2)


t1=1000
tadmix=1001
tsplit=2500
tend=20000
def simple_demography(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture):
    demo=glike.Demo()
    phase1=glike.Phase(0,t1, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1,0],[0,0,1]]), populations=["YRI", "EUR", "HAN"])
    admix=glike.Phase(t1,tadmix, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1-admixture,admixture],[0,0,1]]),populations=["YRI", "EUR", "HAN"])
    phase2=glike.Phase(tadmix,tsplit, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1,0],[0,0,1]]),populations=["YRI", "EUR", "HAN"])
    phase3=glike.Phase(tsplit,tend, [1/NeA,1/NeB],  P=np.array([[1,0], [0,1],[0,0]]),populations=["YRI", "EUR"])
    phase4=glike.Phase(tend,np.inf, [1/NeA],  P=np.array([[1],[0]]), populations=["YRI"])
    demo.add_phase(phase1)
    demo.add_phase(admix)
    demo.add_phase(phase2)
    demo.add_phase(phase3)
    demo.add_phase(phase4)
    
    return demo

demo=simple_demography(t1,tadmix,tsplit, tend, Ne, Ne, Ne, admixture)
demography = miscellaneous.demo_to_demography(demo)
print(demography)

trees=[mod_ts.at(pos).copy() for pos in range(0, 10000000, 10000)]
#trees=[mod_ts.at_index(pos).copy() for pos in range(0, 100000)]#mod_ts.num_trees)]

glike.glike_trees(trees, demo)

def fun(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture):
  demo = simple_demography(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture)
  return glike.glike_trees(trees, demo)

x0 = {"t1":t1,"tadmix":tadmix,"tsplit":tsplit, "tend":tend, "NeA":Ne, "NeB":Ne, "NeC":Ne, "admixture":admixture}
bounds = [(0,tadmix), (t1,tsplit),(tadmix,tend),(tsplit,30000),(0,2*Ne), (0,2*Ne), (0,2*Ne), (0,1)]

#x, logp = estimate.maximize(fun, x0, bounds = bounds)

logl=[]

for l in range(len(trees)):
    dictll={}
               
    dictll.update({"tree":l} )
    dictll.update({"loglikelihood":glike.glike(trees[l], demo)})
    dictll.update({"coal":trees[l].total_branch_length})
    logl.append(dictll)
treell = pd.DataFrame(logl)

plt.plot(treell['tree'], treell['loglikelihood'])
plt.show()

plt.plot(treell['tree'], treell['coal'])
plt.show()


# def find_outliers_IQR(df):

#    q1=df.quantile(0.25)

#    q3=df.quantile(0.75)

#    IQR=q3-q1

#    outliers = df[( (df>(q3+1.5*IQR)))]

#    return outliers

# outliers=find_outliers_IQR(treell['loglikelihood'])

# outTree=trees[outliers.index[0]]

meanll=statistics.mean(treell['loglikelihood'])

def find_adjacent_below_threshold(data, threshold, min_length=2):
    """
    Finds sequences of adjacent values below a threshold in a list.

    Args:
        data (list): The list of values to analyze.
        threshold (float): The threshold value.
        min_length (int): The minimum length of a sequence to be considered.

    Returns:
        list: A list of tuples, each containing the start and end indices of a sequence.
    """

    sequences = []
    start_index = None

    for i, value in enumerate(data):
        if value > threshold:
            if start_index is None:
                start_index = i
        else:
            if start_index is not None and i - start_index >= min_length:
                sequences.append((start_index, i - 1))
            start_index = None

    # Check if the last sequence is still open
    if start_index is not None and len(data) - start_index >= min_length:
        sequences.append((start_index, len(data) - 1))

    return sequences

outliers=find_adjacent_below_threshold(treell['loglikelihood'], meanll, 20)

outtrees=treell[outliers[0][0]:outliers[0][1]]
maxtreeid=outtrees[outtrees['loglikelihood']==max(outtrees['loglikelihood'])]
maxtree=trees[maxtreeid.index[0]]

mcoals=miscellaneous.get_coals(maxtree)

for i in x0:
    print(i)
    print(globals()[i])
    glike.glike(maxtree, demo)
    
    
    

    



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
from IPython.display import SVG

path="/Users/olj5016/Documents/arg_selection/"
#path="/Users/olivia/Documents/arg_selection/"


os.chdir(path)

sys.path.insert(1, "/Users/olj5016/glike/glike/")
#sys.path.insert(1, "/Users/olivia/arg_selection/")
import glike
import estimate
import miscellaneous
import dnisearch
import hsearch


Ne=10000 # unscaled Ne
rr=1e-7 # unscaled recombination rate
mr=1e-8 # unscaled mutation rate
s=0.0 # selection coefficient (0.0 for neutral)
sampleSize=8 # number of indidivuals remebered in tree per generation
selPop=2 # subset value for subpopulation in SLiM 2 for p2, 4 for p22
selTime=17500 # time of selection (generations)
selEnd=20000 # time of selection (generations)
cF = 0.0 # conditional frequency of selected allele (only active when s>0.0)
admixture=0.000000 ## admixture proportion set to 0 to turn admixture off
rep=0 # replicate number

neut = []
for i in range(0,10):
    rep=i
    print(rep)
# set parameter label
    params="{6}_s{0}_sT{1}_sE{2}_sP{3}_cF{4}_admix{5}".format(s,selTime, selEnd,selPop, cF, admixture,rep)

#check if simultion file already ec=xistis, if not simulate demography
    if os.path.isfile("{0}simplegross_{1}.trees".format(path,params))==False:
        # simulate msprime burn in
        burnin = msprime.sim_ancestry(samples=Ne, population_size=Ne, recombination_rate=rr, sequence_length=1e7,   )
        # annotate burn in for slim and export save as .trees file for slim to use to initiate simulation
        burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1, stage="late")
        burnin_ts.dump("{0}burnin_simple_{1}.trees".format(path,params))
        
        #run slim simulation of simple demography with above paramaters
        cmd = "slim -d s=" + str(s) + " -d rep="+str(rep)+" -d admix="+ str(admixture)+" -d sampleSize=" + str(sampleSize)+ " -d selPop=" + str(selPop)+ " -d selTime=" + str(selTime) + " -d cF=" + str(cF)+ " ~/arg_selection/simple_gross.slim"
        print(cmd)
        os.system(cmd)
    
    # load in simulation file
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
    
    # create vector of times indidivuals are saved at
    ind_times = np.unique(ts.individual_times).astype(int)
    
    
    modernInds=ts.samples(time=1) ## sample inds from final generation
    
    mod_ts=ts.simplify(samples=modernInds) # subset tree to just modern individuals

    t1=1000 # t1 (final time)
    tadmix=1001 # time of admixture
    tsplit=2500 # time of split between Eur and Han
    tend=20000 # latest time (start of forward)
    NeA=NeB=NeC=2*Ne # ne of all pops
    
    ## set up glike demography function
    # def simple_demography(tsplit, tend, NeA, NeB, NeC):
    #     demo=glike.Demo()
    #     phase2=glike.Phase(0,tsplit, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1,0],[0,0,1]]))
    #     phase3=glike.Phase(tsplit,tend, [1/NeA,1/NeB],  P=np.array([[1,0], [0,1],[0,1]]))
    #     phase4=glike.Phase(tend,np.inf, [1/NeA],  P=np.array([[1],[1]]))
    #     demo.add_phase(phase2)
    #     demo.add_phase(phase3)
    #     demo.add_phase(phase4)
        
    #     return demo
    def simple_demography(tsplit, tend, NeA, NeB, NeC):

        demo=glike.Demo()

        phase1=glike.Phase(0, tsplit, [1/NeA,1/NeB,1/NeC], populations=["YRI", "EUR", "HAN"])
        phase2=glike.Phase(tsplit,tend, [1/NeA,1/NeB],  P=np.array([[1,0], [0,1],[0,1]]),populations=["YRI", "EUR"])
        phase3=glike.Phase(tend,np.inf, [1/NeA],  P=np.array([[1],[1]]), populations=["YRI"])
    
        demo.add_phase(phase1)
        demo.add_phase(phase2)
        demo.add_phase(phase3)
    
        return demo
    # create demography
    #demo=simple_demography(t1,tadmix,tsplit, tend, Ne, Ne, Ne, admixture)
    demo=simple_demography(tsplit, tend, NeA, NeB, NeC)
    #check demography
    demo.print()
    
    demography = miscellaneous.demo_to_demography(demo)
    print(demography)
    # sample trees from demography
    trees=[mod_ts.at(pos).copy() for pos in range(0, 10000000, 10000)]
    
    ### GLIKE ESTIMATION
    
    # create demography function for maximise function
    # def glike_fun(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture):
    #   demo = simple_demography(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture)
    #   return glike.glike_trees(trees, demo)
    
    
    def glike_fun(tsplit, tend, NeA, NeB, NeC):
      demo = simple_demography(tsplit, tend, NeA, NeB, NeC)
      return glike.glike_trees(trees, demo)
    # list of variables for maximise function
    # glike_x0 = {"t1":t1,"tadmix":tadmix,"tsplit":tsplit, "tend":tend, "NeA":NeA, "NeB":NeB, "NeC":NeC, "admixture":admixture}
    glike_x0 = {"tsplit":tsplit, "tend":tend, "NeA":NeA, "NeB":NeB, "NeC":NeC}

    # bounds for estimate
    # neutral_bounds = [(0,tadmix), (t1,tsplit),(tadmix,tend),(tsplit,np.inf),(0,2*Ne), (0,2*Ne), (0,2*Ne), (0,1)]
    neutral_bounds = [(0,tend),(tsplit,np.inf),(0,10*NeA), (0,10*NeB), (0,10*NeC)]
   
    ## maximise with evenly spaced trees
    est =estimate.maximize(glike_fun, glike_x0, bounds=neutral_bounds)

    ## add estimated parameters to vector to make dt

    neut.append(est[0])
           
# make dt from neutral estimate values          
neutral_est = pd.DataFrame(neut)

#write dt to file
neutral_est.to_string(buf = "neutral_estimates_fixed_1000.txt", index=False)

#calculate averages of neutral parameter estimates
#averages=neutral_est.mean()



import msprime
import pyslim
import os
import sys
import tskit
#import allel
import numpy as np
import pandas as pd
import itertools
# import matplotlib.pyplot as plt
import statistics
from IPython.display import SVG
# import cairosvg

#path to output files
#path="/Users/olj5016/Documents/arg_selection/"
path="/Users/olivia/Documents/arg_selection/"


os.chdir(path)
# path to github repository
#sys.path.insert(1, "/Users/olj5016/glike/glike/")
sys.path.insert(1, "/Users/olivia/arg_selection/")
import glike
import estimate
import miscellaneous
# import dnisearch
# import hsearch


Ne=10000 # unscaled Ne
rr=1e-7 # unscaled recombination rate
mr=1e-8 # unscaled mutation rate
s=0.01 # selection coefficient (0.0 for neutral)
sampleSize=100 # number of indidivuals remebered in tree per generation
selPop=2 # subset value for subpopulation in SLiM 2 for p2, 4 for p22
selTime=19500 # time of selection (generations)
selEnd=20000 # time of selection (generations)
cF = 0.01 # conditional frequency of selected allele (only active when s>0.0)
admixture=0.000000 ## admixture proportion set to 0 to turn admixture off
rep=0 # replicate number

## create parameter label (replicate number_selection coefficient_initation of selection_end of selection_conditional frequency_time to meet conditional frequency_admixture proportion_sample size)

params="{6}_s{0}_sT{1}_sE{2}_sP{3}_cF{4}_admix{5}_sSize{7}".format(s,selTime, selEnd,selPop, cF, admixture,rep,sampleSize)

#check if simultion file already existis, if not simulate demography
if os.path.isfile("{0}simplegross_{1}.trees".format(path,params))==False:
    # simulate msprime burn in
    burnin = msprime.sim_ancestry(samples=Ne, population_size=Ne, recombination_rate=rr, sequence_length=1e7)
    # annotate burn in for slim and export save as .trees file for slim to use to initiate simulation
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1, stage="late")
    burnin_ts.dump("{0}burnin_simple_{1}.trees".format(path,params))

    #run slim simulation of simple demography with above paramaters
    results = "path='" + str(path)+ "'"
    cmd = "slim -d s=" + str(s) + " -d rep="+str(rep)+" -d admix="+ str(admixture)+" -d sampleSize=" + str(sampleSize)+ " -d selPop=" + str(selPop)+ " -d selTime=" + str(selTime) +" -d selEnd=" + str(selEnd)  + " -d cF=" + str(cF)+ ' -d "' + str(results) + '" ~/arg_selection/old_slim.slim' #make sure path to slim file is correct
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
mod_ts=ts.simplify(samples=modernInds) # subset tree to just individuals from selPop


### SWEEP AFTER P2 (C) SPLIT

#overlay mutations on trees
mut_ts = msprime.sim_mutations(mod_ts, rate=mr, discrete_genome=True, keep=True)

def simple_demography(tsplit, tend, NeA, NeB, NeC, tSweep, NeSweep):
      demo=glike.Demo()
      phase1=glike.Phase(0, tSweep, [1/NeA,1/NeB,1/NeSweep], P=np.array([[1,0,0],[0,1,0],[0,0,1]]))
      phase2=glike.Phase(tSweep,tsplit, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1,0],[0,0,1]]))
      phase3=glike.Phase(tsplit,tend, [1/NeA,1/NeB],  P=np.array([[1,0], [0,1],[0,1]]))
      phase4=glike.Phase(tend,np.inf, [1/NeA],  P=np.array([[1],[1]]))
      demo.add_phase(phase1)
      demo.add_phase(phase2)
      demo.add_phase(phase3)
      demo.add_phase(phase4)
      return demo

demo=simple_demography(2501, 20000, 20000, 20000, 20000, 500, 2000)
#set variables
tsplit=2500
tend=20000
NeA=20000
NeB=20000
NeC=20000
NeSweep=20000
tSweep=500


#check demography
demography = miscellaneous.demo_to_demography(demo)
print(demography)
demo.print()

# create population labels for samples
tmp = ["A"] * (2*sampleSize) + ["B"] * (2*sampleSize) + ["C"] * (2*sampleSize) 
samples = {i:pop for i, pop in enumerate(tmp)}

#create gLike function - pass variables and get likelihood of parameter fit
def glike_fun(tsplit, tend, NeA, NeB, NeC, tSweep, NeSweep):
    demo = simple_demography(tsplit, tend, NeA, NeB, NeC, tSweep, NeSweep)
    return glike.glike_trees(trees, demo, samples, kappa=10000)  #return likelihood

glike_x0 = {"tsplit":tsplit, "tend":tend, "NeA":NeA, "NeB":NeB,"NeC":NeC, "tSweep":tSweep, "NeSweep":NeSweep} # dict of initial values to start estimating from

modern_bounds = [(tsplit,tsplit),(tend,tend),(NeA,NeA), (NeB,NeB), (NeC,NeC),(0,tsplit), (0,10*Ne)] # set bounds for estimation values, all fixed expet Ne of interest


win=[] # vec of optimized parameter values
win_est=[] # vec of optimized likelihoods
win_fixed=[] # vec of likelihood under true model
pos=[] # position of window
winsize=50000 #window size
for i in range(0, int(mod_ts.sequence_length),winsize): #loop through the simulated segment
    pos.append(int(i+(winsize/2))) # add position of midpoint of the window to pos vec
    trees=[mut_ts.at(int(i+(winsize/2))).copy()] # extract tree from midpoint of window
    fixed_est= glike.glike_trees(trees, demo, samples,kappa=10000) # calculate likelihood of tree based on fixed true values
    win_fixed.append(fixed_est) # add fixed likelihood estimate to vec
    sel_est=estimate.maximize(glike_fun, glike_x0, bounds=modern_bounds, precision=0.005, epochs=50) #optimise parameters for tree data
    win.append(sel_est[0]) # add optimised param values to vec
    win_est.append(sel_est[1]) # add optimised likelihood to vector

diff=[] # vector of difference between optimised and fixed likelihoods
tdiff=[] # 2 x difference for p-val calc (conducted in R)
for i in range(0,len(win_est)):
    ratio=win_est[i]-win_fixed[i] # calculate differece
    diff.append(ratio) # append diff to vec
    tdiff.append((2*ratio)) # append 2 x diff to vec

# create dt of windowed values
windata = {'Pos': pos, 'l_ratio': tdiff, 'l_fixed': win_fixed, 'l_estimate':win_est}
windows = pd.DataFrame(windata)
# export table
windows.to_csv( "recentwindowedlr_single_{0}.txt".format(params), sep="\t", index=False)

# Run Rscript that will plot tdiff and calculate and plot pvalues for windows
rcmd="Rscript --vanilla /Users/olj5016/arg_selection/pval.R 'recentwindowedlr_single_{0}.txt'".format(params)
os.system(rcmd)

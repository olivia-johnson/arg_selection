import msprime
import pyslim
import os
import tskit
import sys

scaling=10 # scale parameters by 10
Ne=10000 # unscaled Ne
rr=1e-8 # unscaled recombination rate
mr=1e-8 # unscaled mutation rate
s=0.01 # selection coefficient (0.0 for neutral)
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

cmd = "slim -d s=" + str(sS) + " -d sampleSize=" + str(sampleSize)+ " -d selPop=" + str(selPop)+ " -d selTime=" + str(selTime) + " -d cF=" + str(cF)+ " ~/arg_selection/gross_demography.slim"
print(cmd)
os.system(cmd)

ts = tskit.load("/Users/olj5016/Documents/arg_selection/simtree_p{0}_t{1}_s{2}_f{3}_sS{4}.trees".format(selPop, selTime, sS, cF, sampleSize))
mut_ts = msprime.sim_mutations(ts, rate=sM, discrete_genome=False, keep=True)




    
from IPython.display import display
x_limits = [4950000, 5050000]
ts.draw_svg(path="/Users/olj5016/Documents/arg_selection/gross_arg.pdf",y_axis=True, x_lim=x_limits)



import msprime
import pyslim


burnin = msprime.sim_ancestry(samples=1000, population_size=1000, recombination_rate=1e-7, sequence_length=1e7)
burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
burnin_ts.dump("/Users/olj5016/Documents/arg_selection/gross_burnin.trees")



    


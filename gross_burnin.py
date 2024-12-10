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
sampleSize=8 # number of indidivuals remebered in tree per generation
selPop=2 # subset value for subpopulation in SLiM 2 for p2, 4 for p22
selTime=500 # time of selection (generations)
cF = 1 # condntional frequency of selected allele (only active when s>0.0)

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
os.chdir("/Users/olj5016/Documents/arg_selection/GRoSS-master/")
grosscmd="Rscript GRoSS.R -e {0} -d /Users/olj5016/Documents/arg_selection/treeGross.dot -o /Users/olj5016/Documents/arg_selection/gross_out_{1}.tsv".format(filename, params)
os.system(grosscmd)

simpts=mut_ts.simplify()
fts=simpts.first()
fts.draw_svg()
fts.next()
fts.draw_svg()


from datetime import datetime
# names = {"YRI": "African", "CEU": "European", "CHB": "Chinese"}
# colours = {"YRI": "yellow", "CEU": "green", "CHB": "blue"}

# population_map = {p.metadata["id"]: p.id for p in mut_ts.populations()}
sample_populations = list(sorted({simpts.node(u).population for u in simpts.samples()}))
topology_span = {tree.rank(): 0 for tree in tskit.all_trees(len(sample_populations))}

start = datetime.now()
total = 0
for topology_counter, tree in zip(simpts.count_topologies(), simpts.trees()):
    embedded_topologies = topology_counter[sample_populations]
    weight = tree.span / simpts.sequence_length
    for rank, count in embedded_topologies.items():
        topology_span[rank] += count * weight
        total += count
print(f"Counted {total} embedded topologies in {datetime.now() - start} seconds")


# ******************  Simulation Summary  ****************** #

# Specify the output directory for embedded topologies & tree weights
outDIR = '/Users/olj5016/Documents/arg_selection/'
topologies_file_path = os.path.join(outDIR, 'embedded_topologies.pickle')
topologySpan_file_path = os.path.join(outDIR, 'topologySpan.pickle')

population_map = {p.metadata["name"]: p.id for p in mut_ts.populations()}
sample_populations = list(sorted({mut_ts.node(u).population for u in mut_ts.samples()}))
population_map
sample_populations
#' Afr, EEF, Loschbour & Nean
sample_populations = [sample_populations[i] for i in range(len(sample_populations))] # Mbuti, LBK, Losch, Nean


# Identify all tree topologies for sample populations
topology_span = {tree.rank(): 0 for tree in tskit.all_trees(len(sample_populations))}

# Traverse the trees & compute the weights for each tree rank (topolgoy)
from datetime import datetime
start = datetime.now()
total = 0
for topology_counter, tree in zip(mut_ts.count_topologies(), mut_ts.trees()):
    embedded_topologies = topology_counter[sample_populations]
    weight = tree.span / mut_ts.sequence_length
    for rank, count in embedded_topologies.items():
        topology_span[rank] += count * weight
        total += count
print(f"Counted {total} embedded topologies in {datetime.now() - start} seconds")
# >>> Counted 2865712 embedded topologies in 1 day, 1:35:25.138643 seconds

# Save the embedded_topologies collections.Counter object to a file
with open(topologies_file_path, 'wb') as f:
    pickle.dump(embedded_topologies, f)

# Save the topology_span collections.Counter object to a file
with open(topologySpan_file_path, 'wb') as f:
    pickle.dump(topology_span, f)



# Generate Figures of all the topologies & their corresponding weights
names = {"pop_0" : "P0", 
         "0" : "P0", 
         "p1" : "P1",
         "p2" : "P2",
         "p12":"P12",
         "p22":"P22"}
colours = {"pop_0": "black", 
           "p1": "red",
           "p2": "purple",
           "p12":"blue",
           "p22":"orange"}

# Ensure the directory exists
import os
os.makedirs(outDIR, exist_ok=True)
ntips = len(sample_populations)
styles = ".sample text.lab {baseline-shift: super; font-size: 0.7em;}"
node_labels = {}
for p in range(ntips):
    name = mut_ts.population(sample_populations[p]).metadata["name"]
    node_labels[p] = names[name]
    styles += f".n{p}>.sym {{fill: {colours[name]} }}"
total = sum(topology_span.values())
for rank, weight in topology_span.items():
    label = f"{weight/total *100:.1f}% of genome"
    embedded_tree = tskit.Tree.unrank(ntips, rank)
    filename = f"tree_{rank}_{weight}.svg"
    t = embedded_tree.draw_svg(size=(300, 300), style="".join(styles), node_labels=node_labels, x_label=label)
    ## Save the SVG to a file
    #with open(os.path.join(dir_path, filename), 'w') as f:
    #    f.write(t)
    # Save the SVG to a temporary file
    svg_filename = os.path.join(outDIR, "temp.svg")
    with open(svg_filename, 'w') as f:
        f.write(t)

    # Convert SVG to PDF
    pdf_filename = os.path.join(outDIR, filename + ".pdf")
    cairosvg.svg2pdf(url=svg_filename, write_to=pdf_filename)

    # Remove the temporary SVG file
    os.remove(svg_filename)




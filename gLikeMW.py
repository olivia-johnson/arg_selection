gLike Trees 

Williams, Matthew
​
Johnson, Olivia Lenore
​
#********************* >> gLike Analysis
# gLike Kamm Basal Eurasian Model
# msp Kamm Basal Eurasian Simulation

import numpy as np
import tskit
import stdpopsim
import numpy as np
import pandas as pd
import glike
from string import ascii_lowercase
from IPython.display import SVG
import cairosvg
import pickle
import collections
import os

# *** --- *** Read in the tree-sequence
sim_chr = stdpopsim.get_species("HomSap").get_contig("chr21", mutation_rate=stdpopsim.get_species("HomSap").get_demographic_model("AncientEurasia_9K19").mutation_rate).original_coordinates[0]
mts = tskit.load(f"/Users/mkw5910/Documents/aDNA-ARGs/sims/kamm_sims/ts_dir/kammBEur_mts_{sim_chr}.trees")



# *** --- *** Generate a list of trees to feed into gLike

# *** MODE 1: Select speicifc tree indexes 

# Assuming 'mts' is your TreeSequence object
sample_nodes = mts.samples()
for sample_id in sample_nodes:
   node = mts.node(sample_id)
   population_id = node.population
   # Retrieve population metadata
   population = mts.population(population_id)
   population_name = population.metadata.get("name", f"Population {population_id}")
   print(f"Sample {sample_id}: Population {population_id} ({population_name})")

#' Dictionary for the nodes & populations 
population_dict = {}
for sample_id in sample_nodes:
   node = mts.node(sample_id)
   population_id = node.population
   # Retrieve population metadata
   population = mts.population(population_id)
   population_name = population.metadata.get("name", f"Population {population_id}")
   population_dict[sample_id] = population_name
print(population_dict)

# The tree index i want to sample
index_list = [500, 1000, 1500, 2000]
trees = []  # Initialize an empty list
#' Saving the Tree Figures to a file
dir_path = f"/Users/mkw5910/Documents/aDNA-ARGs/sims/kamm_sims/gLike/gLike_kammBEur__msp_KammBEur__mts_{sim_chr}__trees_{len(index_list)}"
os.makedirs(dir_path, exist_ok=True)
# Initialize an empty matrix to store tree interval information
interval_matrix = []
# Start the loop over the tree index list to get the trees and create figures of them.
for index in index_list:
   try:
       tree_at_index = mts.at_index(index)
       # *** --------------- Generate PDF for each tree
       t = tree_at_index.draw_svg(size=(2000, 500), node_labels=population_dict, y_axis=True)
       filename = f"msp_kammBEur__index_{index}.svg"
       ## Save the SVG to a file
       # Save the SVG to a temporary file
       svg_filename = os.path.join(dir_path, "temp.svg")
       with open(svg_filename, 'w') as f:
           f.write(t)
       # Convert SVG to PDF
       pdf_filename = os.path.join(dir_path, filename + ".pdf")
       cairosvg.svg2pdf(url=svg_filename, write_to=pdf_filename)
       # Remove the temporary SVG file
       os.remove(svg_filename)

       # *** --------------- Continue with creating the list of trees
       trees.append(tree_at_index)  # Append the tree to the list
       # Extract interval information
       interval_start = tree_at_index.interval.left
       interval_end = tree_at_index.interval.right
       # Append interval data to the matrix
       interval_matrix.append([tree_at_index.interval, interval_start, interval_end])
       # Print interval details
       print(f"Tree {index} covers interval: {tree_at_index.interval} "
             f"(Length: {interval_end - interval_start:.2f} units)")
   except IndexError:
       print(f"Index {index} is out of range for the given TreeSequence.")

# Now 'interval_matrix' contains the desired information
print(f"Number of trees in the list: {len(trees)}")
print(interval_matrix)




# *** MODE 2: Sample X trees between positions Y and Z spaced by length L
seqL = mts.sequence_length
list_of_trees = [mts.at(pos).copy() for pos in range(1, int(seqL), 500000)] 
print(f"Number of trees in the list: {len(list_of_trees)}")
#list_of_trees = [mts.at(pos).copy() for pos in range(3000000, 30000000, 1350000)] # 20 trees
#print(f"Number of trees in the list: {len(list_of_trees)}")



# gLike Demographic Model 
def kamm_BEur(N_Mbuti, N_WHG_Basal_MA1_Ust, 
                  N_EEF, N_Sard, N_Losch, 
                  N_Han, N_Nean_a, 
                  N_WHG_Mbuti, N_WHG_Han, 
                  N_ANC, 

                  p_WHG_Sard, p_Basal_EEF, p_Nean_Eur,

                  Tadmix_WHG_Sard, Tsplit_EEF_Sard, Tadmix_Basal_EEF,
                  Tsplit_WHG_EEF, Tsplit_WHG_MA1, Tsplit_WHG_Han,
                  Tsplit_WHG_Ust, Tadmix_Nean_Eur, Tsplit_WHG_Basal,
                  Tsplit_WHG_Mbuti, Tsplit_WHG_Nean):
 demo = glike.Demo()
 # Time:       0: Tadmix_WHG_Sard = (1230 gens)
 # Start:      N=9 (1=Mbuti, 2=BasalE, 3=EEF, 4=Sard, 5=Losch, 6=MA1, 7=Han, 8=Ust, 9=Nean)
 # End:        N=9 (1=Mbuti, 2=BasalE, 3=EEF, 4=Sard, 5=Losch, 6=MA1, 7=Han, 8=Ust, 9=Nean)
 demo.add_phase(glike.Phase(0, Tadmix_WHG_Sard, 
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_EEF, 1/N_Sard, 
                             1/N_Losch, 1/N_WHG_Basal_MA1_Ust, 1/N_Han, 1/N_WHG_Basal_MA1_Ust, 1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0, 0, 0, 0, 0, 0], 
                                          [0, 1, 0, 0, 0, 0, 0, 0, 0], 
                                          [0, 0, 1, 0, 0, 0, 0, 0, 0],
                                          [0, 0, 0, 1, 0, 0, 0, 0, 0], 
                                          [0, 0, 0, 0, 1, 0, 0, 0, 0],
                                          [0, 0, 0, 0, 0, 1, 0, 0, 0],
                                          [0, 0, 0, 0, 0, 0, 1, 0 , 0],
                                          [0, 0, 0, 0, 0, 0, 0, 1, 0],
                                          [0, 0, 0, 0, 0, 0, 0, 0, 1]])))
 # Time:       Tadmix_WHG_Sard = (1230 gens) : Tsplit_EEF_Sard = (7690 gens)
 # Start:      N=9 (1=Mbuti, 2=BasalE, 3=EEF, 4=Sard, 5=Losch, 6=MA1, 7=Han, 8=Ust, 9=Nean)
 # End:        N=9 (1=Mbuti, 2=BasalE, 3=EEF, 4=Sard, 5=Losch, 6=MA1, 7=Han, 8=Ust, 9=Nean)
 # Migration:  Losch -> Sard (p_WHG_Sard = 0.0317)
 demo.add_phase(glike.Phase(Tadmix_WHG_Sard, Tsplit_EEF_Sard,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_EEF, 1/N_Sard, 
                             1/N_Losch, 1/N_WHG_Basal_MA1_Ust, 1/N_Han, 1/N_WHG_Basal_MA1_Ust, 1/N_Nean_a],
                            P = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0], 
                                          [0, 1, 0, 0, 0, 0, 0, 0, 0], 
                                          [0, 0, 1, 0, 0, 0, 0, 0, 0],
                                          [0, 0, 0, 1-p_WHG_Sard, p_WHG_Sard, 0, 0, 0, 0], # Losch (p)-> Sard (forwards)
                                          [0, 0, 0, 0, 1, 0, 0, 0, 0],
                                          [0, 0, 0, 0, 0, 1, 0, 0, 0],
                                          [0, 0, 0, 0, 0, 0, 1, 0 , 0],
                                          [0, 0, 0, 0, 0, 0, 0, 1, 0],
                                          [0, 0, 0, 0, 0, 0, 0, 0, 1]])))
 # Time:       Tsplit_EEF_Sard = (7690 gens) : Tadmix_Basal_EEF = (33700 gens)
 # Start:      N=9 (1=Mbuti, 2=BasalE, 3=EEF, 4=Sard, 5=Losch, 6=MA1, 7=Han, 8=Ust, 9=Nean)
 # End:        N=8 (1=Mbuti, 2=BasalE, 3=EEF, 4=Losch, 5=MA1, 6=Han, 7=Ust, 8=Nean)
 # Migration:  Sard -> LBK (p=1) 
 demo.add_phase(glike.Phase(Tsplit_EEF_Sard, Tadmix_Basal_EEF,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_EEF, 
                             1/N_Losch, 1/N_WHG_Basal_MA1_Ust, 1/N_Han, 1/N_WHG_Basal_MA1_Ust, 1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0, 0, 0, 0, 0], 
                                          [0, 1, 0, 0, 0, 0, 0, 0], 
                                          [0, 0, 1, 0, 0, 0, 0, 0],
                                          [0, 0, 1, 0, 0, 0, 0, 0], # Sard split from EEF
                                          [0, 0, 0, 1, 0, 0, 0, 0], 
                                          [0, 0, 0, 0, 1, 0, 0, 0],
                                          [0, 0, 0, 0, 0, 1, 0, 0],
                                          [0, 0, 0, 0, 0, 0, 1, 0 ],
                                          [0, 0, 0, 0, 0, 0, 0, 1]])))
 # Time:       Tadmix_Basal_EEF = (33700 gens) : Tsplit_WHG_EEF = (37700 gens)
 # Start:      N=8 (1=Mbuti, 2=BasalE, 3=EEF, 4=Losch, 5=MA1, 6=Han, 7=Ust, 8=Nean)
 # End:        N=8 (1=Mbuti, 2=BasalE, 3=EEF, 4=Losch, 5=MA1, 6=Han, 7=Ust, 8=Nean)
 # Migration:  BasalE -> LBK (p_Basal_EEF = 0.0936)
 demo.add_phase(glike.Phase(Tadmix_Basal_EEF, Tsplit_WHG_EEF,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_EEF, 
                             1/N_Losch, 1/N_WHG_Basal_MA1_Ust, 1/N_Han, 1/N_WHG_Basal_MA1_Ust, 1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0, 0, 0, 0, 0], 
                                          [0, 1, 0, 0, 0, 0, 0, 0], 
                                          [0, p_Basal_EEF, 1-p_Basal_EEF, 0, 0, 0, 0, 0], # BasalE (p)-> EEF (forwards)
                                          [0, 0, 0, 1, 0, 0, 0, 0], 
                                          [0, 0, 0, 0, 1, 0, 0, 0],
                                          [0, 0, 0, 0, 0, 1, 0, 0],
                                          [0, 0, 0, 0, 0, 0, 1, 0],
                                          [0, 0, 0, 0, 0, 0, 0, 1]])))
 # Time:       Tsplit_WHG_EEF = (37700 gens) : Tsplit_WHG_MA1 = (44900 gens)
 # Start:      N=8 (1=Mbuti, 2=BasalE, 3=EEF, 4=Losch, 5=MA1, 6=Han, 7=Ust, 8=Nean)
 # End:        N=7 (1=Mbuti, 2=BasalE, 3=Losch, 4=MA1, 5=Han, 6=Ust, 7=Nean)
 # Migration:  EEF -> WHG (p=1)
 demo.add_phase(glike.Phase(Tsplit_WHG_EEF, Tsplit_WHG_MA1,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_Losch, 
                             1/N_Han, 1/N_WHG_Basal_MA1_Ust, 1/N_WHG_Basal_MA1_Ust, 1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0, 0, 0, 0], 
                                          [0, 1, 0, 0, 0, 0, 0], 
                                          [0, 0, 1, 0, 0, 0, 0], # EEF split from Losch
                                          [0, 0, 1, 0, 0, 0, 0], 
                                          [0, 0, 0, 1, 0, 0, 0],
                                          [0, 0, 0, 0, 1, 0, 0],
                                          [0, 0, 0, 0, 0, 1, 0],
                                          [0, 0, 0, 0, 0, 0, 1]])))
 # Time:       Tsplit_WHG_MA1 = (44900 gens) : Tsplit_WHG_Han = (50400 gens)
 # Start:      N=7 (1=Mbuti, 2=BasalE, 3=Losch, 4=MA1, 5=Han, 6=Ust, 7=Nean)
 # End:        N=6 (1=Mbuti, 2=BasalE, 3=Losch, 4=Han, 5=Ust, 6=Nean)
 # Migration:  MA1 -> WHG (p=1)
 demo.add_phase(glike.Phase(Tsplit_WHG_MA1, Tsplit_WHG_Han,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_Losch, 
                             1/N_Han, 1/N_WHG_Basal_MA1_Ust, 1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0, 0, 0], 
                                          [0, 1, 0, 0, 0, 0], 
                                          [0, 0, 1, 0, 0, 0],
                                          [0, 0, 1, 0, 0, 0], # MA1 split from Losch
                                          [0, 0, 0, 1, 0, 0],
                                          [0, 0, 0, 0, 1, 0],
                                          [0, 0, 0, 0, 0, 1]])))
 # Time:       Tsplit_WHG_Han = (50400 gens) : Tsplit_WHG_Ust = (51500 gens)
 # Start:      N=6 (1=Mbuti, 2=BasalE, 3=Losch, 4=Han, 5=Ust, 6=Nean)
 # End:        N=5 (1=Mbuti, 2=BasalE, 3=Losch, 4=Ust, 5=Nean)
 # Migration:  Han -> WHG (p=1)
 demo.add_phase(glike.Phase(Tsplit_WHG_Han, Tsplit_WHG_Ust,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_WHG_Han, 
                             1/N_WHG_Basal_MA1_Ust, 1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0, 0], 
                                          [0, 1, 0, 0, 0], 
                                          [0, 0, 1, 0, 0],
                                          [0, 0, 1, 0, 0], # Han split from Losch
                                          [0, 0, 0, 1, 0],
                                          [0, 0, 0, 0, 1]])))
 # Time:       Tsplit_WHG_Ust = (51500 gens) : Tadmix_Nean_Eur = (56800 gens)
 # Start:      N=5 (1=Mbuti, 2=BasalE, 3=Losch, 4=Ust, 5=Nean)
 # End:        N=4 (1=Mbuti, 2=BasalE, 3=Losch, 4=Nean)
 # Migration:  Ust -> WHG (p=1)
 demo.add_phase(glike.Phase(Tsplit_WHG_Ust, Tadmix_Nean_Eur,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_WHG_Han, 
                              1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0], 
                                          [0, 1, 0, 0], 
                                          [0, 0, 1, 0],
                                          [0, 0, 1, 0], # Ust split from Losch
                                          [0, 0, 0, 1]])))
 # Time:       Tadmix_Nean_Eur = (56800 gens) : Tsplit_WHG_Basal = (79800 gens)
 # Start:      N=4 (1=Mbuti, 2=BasalE, 3=Losch, 4=Nean)
 # End:        N=4 (1=Mbuti, 2=BasalE, 3=Losch, 4=Nean)
 # Migration:  Nean -> Eur (p_Nean_Eur = 0.0296)
 demo.add_phase(glike.Phase(Tadmix_Nean_Eur, Tsplit_WHG_Basal,
                            [1/N_Mbuti, 1/N_WHG_Basal_MA1_Ust, 1/N_WHG_Han, 
                              1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0, 0], 
                                          [0, 1, 0, 0], 
                                          [0, 0, 1-p_Nean_Eur, p_Nean_Eur], # Nean (p)-> Losch (forward)
                                          [0, 0, 0, 1]])))
 # Time:       Tsplit_WHG_Basal = (79800 gens) : Tsplit_WHG_Mbuti = (95800 gens)
 # Start:      N=4 (1=Mbuti, 2=BasalE, 3=Losch, 4=Nean)
 # End:        N=3 (1=Mbuti, 2=Losch, 3=Nean)
 # Migration:  BasalE -> WHG (p=1)
 demo.add_phase(glike.Phase(Tsplit_WHG_Basal, Tsplit_WHG_Mbuti,
                            [1/N_Mbuti, 1/N_WHG_Han, 
                              1/N_Nean_a],
                             P = np.array([
                                          [1, 0, 0], 
                                          [0, 1, 0], # Basal split from Losch
                                          [0, 1, 0],
                                          [0, 0, 1]])))
 # Time:       Tsplit_WHG_Mbuti = (95800 gens) : Tsplit_WHG_Nean = (696000 gens)
 # Start:      N=3 (1=Mbuti, 2=Losch, 3=Nean)
 # End:        N=2 (1=Losch, 2=Nean)
 # Migration:  WHG -> Mbuti (p=1)
 demo.add_phase(glike.Phase(Tsplit_WHG_Mbuti, Tsplit_WHG_Nean,
                            [1/N_WHG_Mbuti,
                              1/N_Nean_a],
                             P = np.array([
                                          [1, 0], # Mbuti split from Losch 
                                          [1, 0],
                                          [0, 1]])))
 # Time:       Tsplit_WHG_Nean = (696000 gens) : Inf
 # Start:      N=2 (1=Mbuti, 2=Nean)
 # End:        N=1 (1=Mbuti)
 # Migration:  Nean -> Mbuti (p=1)
 demo.add_phase(glike.Phase(Tsplit_WHG_Nean, np.inf,
                            [N_ANC],
                            P = np.array([[1], 
                                          [1]]))) # Nean split from Losch
 return demo

                            
                           
# Estimating parameters
names = ["N_Mbuti", "N_EEF", "N_Sard", "N_WHG_Basal_MA1_Ust", "N_Han", "N_Nean_a",
        "N_WHG_Han", "N_WHG_Mbuti", "N_ANC", "N_Losch", 

        "p_WHG_Sard", "p_Basal_EEF", "p_Nean_Eur",

        "Tadmix_WHG_Sard",  "Tsplit_EEF_Sard", "Tadmix_Basal_EEF", 
        "Tsplit_WHG_EEF", "Tsplit_WHG_MA1", "Tsplit_WHG_Han", "Tsplit_WHG_Ust", "Tadmix_Nean_Eur",
        "Tsplit_WHG_Basal", "Tsplit_WHG_Mbuti", "Tsplit_WHG_Nean"]
values = [17300, 75.7, 15000, 1920, 6300, 86.9, 2340, 29100, 18200, 1200, 
         
         0.0317, 0.0936,0.0296,

         1230, 7690, 33700,  37700, 44900, 50400, 51500, 
         56800, 79800, 95800, 696000]

limits = [(100,100000), (100,100000), (100,100000), (100,100000), (100,100000), 
         (100,100000), (100,100000), (100,100000), (100,100000), (100,100000), 
         
         (0.001,0.999), (0.001,0.999), (0.001,0.999), 

         (5, "Tadmix_WHG_Sard"), 
         ("Tsplit_EEF_Sard", "Tadmix_Basal_EEF"),
         ("Tsplit_EEF_Sard", "Tsplit_WHG_EEF"),
         ("Tadmix_Basal_EEF", "Tsplit_WHG_MA1"),
         ("Tsplit_WHG_EEF", "Tsplit_WHG_Han"),
         ("Tsplit_WHG_MA1", "Tsplit_WHG_Ust"),
         ("Tsplit_WHG_Han", "Tadmix_Nean_Eur"),
         ("Tsplit_WHG_Ust", "Tsplit_WHG_Basal"),
         ("Tadmix_Nean_Eur", "Tsplit_WHG_Mbuti"),
         ("Tsplit_WHG_Basal", "Tsplit_WHG_Nean"),
         ("Tsplit_WHG_Mbuti", 10e5)]

search = glike.Search(names, values, limits, precision = 0.02)

x, logp = glike.estimate(list_of_trees, kamm_BEur, search, prune = 0.5)


# *** --- *** Write Results to File
# Parameter names
parameters = ['N_Mbuti', 'N_WHG_Basal_MA1_Ust', 'N_EEF', 'N_Sard', 'N_Losch', 'N_Han', 'N_Nean_a', 'N_WHG_Mbuti', 'N_WHG_Han', 'N_ANC', 
             'p_WHG_Sard', 'p_Basal_EEF', 'p_Nean_Eur', 'Tadmix_WHG_Sard', 'Tsplit_EEF_Sard', 'Tadmix_Basal_EEF','Tsplit_WHG_EEF', 'Tsplit_WHG_MA1', 'Tsplit_WHG_Han',
             'Tsplit_WHG_Ust', 'Tadmix_Nean_Eur', 'Tsplit_WHG_Basal', 'Tsplit_WHG_Mbuti', 'Tsplit_WHG_Nean']

# True values
values = [17300, 75.7, 15000, 1920, 6300, 86.9, 2340, 29100, 18200, 1200, 
         0.0317, 0.0936,0.0296,
         1230, 7690, 33700,  37700, 44900, 50400, 51500, 
         56800, 79800, 95800, 696000]
# Estimated values
x 
# Model likelihood
logp

# Create a Summary DataFrame
df = pd.DataFrame([values, x], columns=parameters, index=['Simulated value', 'Estimated value'])
# Subtract the 'Simulated value' row from the 'Estimated value' row
df.loc['Estimated - Simulated'] = df.loc['Estimated value'] - df.loc['Simulated value'] 
# Add a new column 'Model Likelihood' with the value of logp in the second row
df['Model Likelihood'] = [None, logp, None]
# Your updated dataframe will look like this:
print(df)


# Output directory and filename
outDIR = f"/Users/mkw5910/Documents/aDNA-ARGs/sims/kamm_sims/gLike/params_gLike_kammBEur__msp_KammBEur__mts_{sim_chr}__trees_{len(list_of_trees)}.txt"

# Write the DataFrame to a tab-separated file
df.to_csv(outDIR, sep='\t')






#' Kamm Model Parameter Variables
N_Mbuti = 17300
N_EEF = 75.7
N_Sard = 15000
N_WHG_Basal_MA1_Ust = 1920
N_Han = 6300
N_Nean_a = 86.9
N_WHG_Han = 2340
N_WHG_Mbuti = 29100
N_ANC = 18200
N_Losch = 1200

p_WHG_Sard = 0.0317
p_Basal_EEF = 0.0936
p_Nean_Eur = 0.0296

Tadmix_WHG_Sard = 1230
Tsplit_EEF_Sard = 7690
Tadmix_Basal_EEF = 33700
Tsplit_WHG_EEF = 37700
Tsplit_WHG_MA1 = 44900
Tsplit_WHG_Han = 50400
Tsplit_WHG_Ust = 51500
Tadmix_Nean_Eur = 56800
Tsplit_WHG_Basal = 79800
Tsplit_WHG_Mbuti = 95800
Tsplit_WHG_Nean = 696000

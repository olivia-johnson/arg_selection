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
import cairosvg

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
s=1.0 # selection coefficient (0.0 for neutral)
sampleSize=100 # number of indidivuals remebered in tree per generation
selPop=1 # subset value for subpopulation in SLiM 2 for p2, 4 for p22
selTime=17000 # time of selection (generations)
selEnd=20000 # time of selection (generations)
cF = 1.0 # conditional frequency of selected allele (only active when s>0.0)
cFTime = 17499
admixture=0.000000 ## admixture proportion set to 0 to turn admixture off
rep=0 # replicate number

sel = []
for i in range(0,10):
    rep=i
    print(rep)
# set parameter label
    params="{6}_s{0}_sT{1}_sE{2}_sP{3}_cF{4}_cFT{8}_admix{5}_sSize{7}".format(s,selTime, selEnd,selPop, cF, admixture,rep,sampleSize,cFTime)
    
    #check if simultion file already ec=xistis, if not simulate demography
    if os.path.isfile("{0}simplegross_{1}.trees".format(path,params))==False:
        # simulate msprime burn in
        burnin = msprime.sim_ancestry(samples=Ne, population_size=Ne, recombination_rate=rr, sequence_length=1e7)
        # annotate burn in for slim and export save as .trees file for slim to use to initiate simulation
        burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1, stage="late")
        burnin_ts.dump("{0}burnin_simple_{1}.trees".format(path,params))
        
        #run slim simulation of simple demography with above paramaters
        cmd = "slim -d s=" + str(s) + " -d rep="+str(rep)+" -d admix="+ str(admixture)+" -d sampleSize=" + str(sampleSize)+ " -d selPop=" + str(selPop)+ " -d selTime=" + str(selTime) +" -d selEnd=" + str(selEnd) + " -d cF=" + str(cF)+" -d cFTime=" + str(cFTime)+ " ~/arg_selection/simple_gross.slim"
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
    mp2=ts.samples(selPop,time=1) #sample inds from selpop and final generation
    sel1=np.append(mp2, (260,261)) # individuals from selected population and one individual from a neutral population
    mod_ts=ts.simplify(samples=modernInds) # subset tree to just indidivuals from selPop
    
    #overlay mutations on trees
    mut_ts = msprime.sim_mutations(mod_ts, rate=mr, discrete_genome=True, keep=True)
    
    # import demography values from estimates of neutral trees
    # neutral_est = pd.read_csv("neutral_estimates_recalc_splitmigration.txt", sep="\t")
    # # calculate averages of the neutral values to use in demograpy
    # averages=neutral_est.mean()
    
    
    # ## Assign neutral averages to dmeography parameters 
    # # t1=int(averages["t1"]) # t1 (final time)
    # # tadmix=int(averages["tadmix"]) # time of admixture
    # tsplit=int(averages["tsplit"]) # time of split between Eur and Han
    # tend=int(averages["tend"]) # latest time (start of forward when pop A split from the ancestory of B and C)
    # NeA=int(averages["NeA"]) # population A (p0) Ne
    # NeB=int(averages["NeB"])# population B (p1) Ne
    # NeC=int(averages["NeC"])# population C (p2) Ne
    # # admixture=averages["admixture"]  #admixture proportion
    
    ## set up glike demography function
    # def simple_demography(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture):
    #     demo=glike.Demo()
    #     phase1=glike.Phase(0,t1, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1,0],[0,0,1]]), populations=["YRI", "EUR", "HAN"])
    #     admix=glike.Phase(t1,tadmix, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1-admixture,admixture],[0,0,1]]),populations=["YRI", "EUR", "HAN"])
    #     phase2=glike.Phase(tadmix,tsplit, [1/NeA,1/NeB,1/NeC], P=np.array([[1,0,0],[0,1,0],[0,0,1]]),populations=["YRI", "EUR", "HAN"])
    #     phase3=glike.Phase(tsplit,tend, [1/NeA,1/NeB],  P=np.array([[1,0], [0,1],[0,0]]),populations=["YRI", "EUR"])
    #     phase4=glike.Phase(tend,np.inf, [1/NeA],  P=np.array([[1],[0]]), populations=["YRI"])
    #     demo.add_phase(phase1)
    #     demo.add_phase(admix)
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
    # demo=simple_demography(t1,tadmix,tsplit, tend, Ne, Ne, Ne, admixture)
    # demo=simple_demography(tsplit, tend, NeA, NeB, NeC)
    
    tsplit=2500
    tend=20000
    NeA=20000
    NeB=20000
    NeC=20000
    demo=simple_demography(2500, 20000, 20000, 20000, 20000)
    
    
    #check demography
    demography = miscellaneous.demo_to_demography(demo)
    print(demography)
    # sample trees from demography
    trees=[mut_ts.at(pos).copy() for pos in range(0, 10000000, 100000)]
    
    #generate overall likelihood value (will bea meanll*10)
    glike.glike_trees(trees, demo) 
    
    
    # # plot trees
    # index_list = [20, 5000000]
    # plottrees = []  # Initialize an empty list
    # #' Saving the Tree Figures to a file
    # dir_path = f"/Users/olj5016/Documents/arg_selection/gliketrees"
    # os.makedirs(dir_path, exist_ok=True)
    # # Initialize an empty matrix to store tree interval information
    # interval_matrix = []
    
    # styles = []
    # # Create a style for each population, programmatically (or just type the string by hand)
    # for colour, p in zip(['red', 'green', 'blue'], mod_ts.populations()):
    #     # target the symbols only (class "sym")
    #     st = f".node.p{p.id} > .sym " + "{" + f"fill: {colour}" + "}"
    #     styles.append(st)
    #     print(f'"{st}" applies to nodes from population {p.metadata["name"]} (id {p.id})')
    # css_string = " ".join(styles)
    # print(f'CSS string applied:\n    "{css_string}"')
    # # Start the loop over the tree index list to get the trees and create figures of them.
    # for index in index_list:
    #        tree_at_index = mod_ts.at(index)
    #        # *** --------------- Generate PDF for each tree
    #        t = tree_at_index.draw_svg(size=(2000, 500), node_labels={},    # Remove all node labels for a clearer viz
    #        style=css_string, y_axis=True)
    #        filename = f"glike_{params}_all_{index}.svg"
    #        ## Save the SVG to a file
    #        # Save the SVG to a temporary file
    #        svg_filename = os.path.join(dir_path, "temp.svg")
    #        with open(svg_filename, 'w') as f:
    #            f.write(t)
    #        # Convert SVG to PDF
    #        pdf_filename = os.path.join(dir_path, filename + ".pdf")
    #        cairosvg.svg2pdf(url=svg_filename, write_to=pdf_filename)
    #        # Remove the temporary SVG file
    #        os.remove(svg_filename)
    
    ### GLIKE ESTIMATION
    
    # create demography function for maximise function
    # def glike_fun(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture):
    #   demo = simple_demography(t1,tadmix,tsplit, tend, NeA, NeB, NeC, admixture)
    #   return glike.glike_trees(trees, demo)
    
    tmp = ["YRI"] * (2*sampleSize) + ["EUR"] * (2*sampleSize) + ["HAN"] * (2*sampleSize) 
    samples = {i:pop for i, pop in enumerate(tmp)}
    tmpp2 = ["HAN"] * (2*sampleSize) 
    samplesp2 = {i:pop for i, pop in enumerate(tmpp2)}
    
    
    def glike_fun(tsplit, tend, NeA, NeB, NeC):
       demo = simple_demography(tsplit, tend, NeA, NeB, NeC)
       return glike.glike_trees(trees, demo, samples, kappa=10000)
    
    
    # list of variables for maximise function
    # glike_x0 = {"t1":t1,"tadmix":tadmix,"tsplit":tsplit, "tend":tend, "NeA":NeA, "NeB":NeB, "NeC":NeC, "admixture":admixture}
    # glike_x0 = {"tsplit":2000, "tend":15000, "NeA":15000, "NeB":15000, "NeC":15000}
    # glike_x0 = {"tsplit":tsplit, "tend":tend, "NeA":NeA, "NeB":NeB, "NeC":NeC}
    glike_x0 = {"tsplit":tsplit, "tend":tend, "NeA":NeA, "NeB":NeB, "NeC":20000}

    # list of bound for maximise function (only works if set for the selPop)
    # glike_bounds = [(t1,t1), (tadmix,tadmix),(tsplit,tsplit),(tend,tend),(NeA,NeA), (NeB,NeB), (0,10*Ne), (admixture,admixture)]
    glike_bounds = [(tsplit,tsplit),(tend,tend),(NeA,NeA), (NeB,NeB), (0,10*Ne)]
    neutral_bounds = [(0,tend),(tsplit,np.inf),(0,10*NeA), (0,10*NeB), (0,10*NeC)]
    
    # calculate branch length and loglikelihood for each simulated tree
#     logl=[]
    
#     for l in range(len(trees)):
#         dictll={}
                   
#         dictll.update({"tree":l} )
#         dictll.update({"loglikelihood":glike.glike(trees[l], demo,samples)})
#         dictll.update({"branch_length":trees[l].total_branch_length})
#         logl.append(dictll)
#     treell = pd.DataFrame(logl)
    
#     ## plot likelihood and branch length
#     plt.plot(treell['tree'], treell['loglikelihood'])
#     plt.savefig('likelihood_p2_{0}.pdf'.format(params))
#     plt.show()
    
#     plt.plot(treell['tree'], treell['branch_length'])
#     plt.savefig('branchlength_p2_{0}.pdf'.format(params))
#     plt.show()
    
#     ## maximise with evenly spaced trees
#     sel_est=estimate.maximize(glike_fun, glike_x0, bounds=glike_bounds)
#     sel.append(sel_est[0])
# sel_est = pd.DataFrame(sel)




##### Windowed glike  #####

win=[]
pos=[]
for i in range(0, int(mod_ts.sequence_length),10000):
    pos.append(i)
    trees=[mut_ts.at(pos).copy() for pos in range(i, (i+20))]
    sel_est=estimate.maximize(glike_fun, glike_x0, bounds=glike_bounds)
    win.append(sel_est[0])

windows = pd.DataFrame(win)
windows['Pos']=pos

plt.plot(windows['Pos'], windows['NeC'])
plt.savefig('windowed_{0}.pdf'.format(params))
plt.show()
## Multiple treees per window
win=[]
win_est=[]
win_fixed=[]
pos=[]
winsize=50000
for i in range(0, int(mod_ts.sequence_length),winsize):
    pos.append(i)
    trees=[mut_ts.at(pos).copy() for pos in range(i, (i+winsize), int(winsize/20))]
    fixed_est= glike.glike_trees(trees, demo, samples, kappa=10000)
    win_fixed.append(fixed_est)
    sel_est=estimate.maximize(glike_fun, glike_x0, bounds=glike_bounds, precision =0.005, epochs=50)
    win.append(sel_est[0])
    win_est.append(sel_est[1])

diff=[]
tdiff=[]
for i in range(0,len(win_est)):
    ratio=win_est[i]-win_fixed[i]
    diff.append(ratio)
    tdiff.append((2*ratio))
    
windata = {'Pos': pos, 'l_ratio': tdiff, 'l_fixed': win_fixed, 'l_estimate':win_est}
windows = pd.DataFrame(windata)

windows.to_csv( "windowedlr_mulitple_s{0}_sT{1}_sE{2}_sP{3}_cF{4}_admix{5}.txt".format(s,selTime, selEnd,selPop, cF, admixture), sep="\t", index=False)

    
# plt.plot(range(0,len(win_est)), win_est)

# plt.plot(pos, tdiff)
# plt.savefig('2diff_multiple_k10000_{0}.pdf'.format(params))
# plt.show

# plt.plot(pos, win_est)
# plt.savefig('estlikeli_multiple_k10000_{0}.pdf'.format(params))
# plt.show


# plt.plot(pos, win_fixed)
# plt.savefig('fixlikeli_multiple_k10000_{0}.pdf'.format(params))
# plt.show


## Single tree per window
win=[]
win_est=[]
win_fixed=[]
pos=[]
winsize=50000
for i in range(0, int(mod_ts.sequence_length),winsize):
    pos.append(int(i+(winsize/2)))
    trees=[mut_ts.at(int(i+(winsize/2))).copy()]
    fixed_est= glike.glike_trees(trees, demo, samples,kappa=10000)
    win_fixed.append(fixed_est)
    sel_est=estimate.maximize(glike_fun, glike_x0, bounds=glike_bounds, precision=0.005, epochs=50)
    win.append(sel_est[0])
    win_est.append(sel_est[1])

diff=[]
tdiff=[]
for i in range(0,len(win_est)):
    ratio=win_est[i]-win_fixed[i]
    diff.append(ratio)
    tdiff.append((2*ratio))


windata = {'Pos': pos, 'l_ratio': tdiff, 'l_fixed': win_fixed, 'l_estimate':win_est}
windows = pd.DataFrame(windata)

windows.to_csv( "windowedlr_single_{0}.txt".format(params), sep="\t", index=False)

# plt.plot(pos, tdiff)
# plt.savefig('2diff_singletree_k10000_{0}.pdf'.format(params))
# plt.show


# plt.plot(pos, win_est)
# plt.savefig('estlikeli_singletree_k10000_{0}.pdf'.format(params))
# plt.show


# plt.plot(pos, win_fixed)
# plt.savefig('fixlikeli_singletree_k10000_{0}.pdf'.format(params))
# plt.show


# ### GLIKE ESTIMATES FOR OUTLIER TREES ###

# # calculate mean likelihood of trees
# meanll=statistics.mean(treell['loglikelihood'])


# #function to find regions with trees over the threshold (meanll)
# def find_adjacent_above_threshold(data, threshold, min_length=2):
#     """
#     Finds sequences of adjacent values above a threshold in a list.

#     Args:
#         data (list): The list of values to analyze.
#         threshold (float): The threshold value.
#         min_length (int): The minimum length of a sequence to be considered.

#     Returns:
#         list: A list of tuples, each containing the start and end indices of a sequence.
#     """

#     sequences = []
#     start_index = None

#     for i, value in enumerate(data):
#         if value > threshold:
#             if start_index is None:
#                 start_index = i
#         else:
#             if start_index is not None and i - start_index >= min_length:
#                 sequences.append((start_index, i - 1))
#             start_index = None

#     # Check if the last sequence is still open
#     if start_index is not None and len(data) - start_index >= min_length:
#         sequences.append((start_index, len(data) - 1))

#     return sequences


# ## identify trees over the threshold
# outliers=find_adjacent_above_threshold(treell['loglikelihood'], meanll,30)

# ## extract outliers into new trees vector
# outtrees=[]
# for i in range(len(outliers)):
#     outtrees=outtrees+list(range(outliers[i][0],outliers[i][1]))
# treell
# trees=[mod_ts.at(pos*1000).copy() for pos in outtrees]

# tlogl=[]

# for l in range(len(trees)):
#     tdictll={}
               
#     tdictll.update({"tree":l} )
#     tdictll.update({"loglikelihood":glike.glike(trees[l], demo, samples)})
#     tdictll.update({"branch_length":trees[l].total_branch_length})
#     tlogl.append(tdictll)
# ttreell = pd.DataFrame(tlogl)

# ## plot likelihood and branch length
# plt.plot(ttreell['tree'], ttreell['loglikelihood'])
# plt.show()

# plt.plot(ttreell['tree'], ttreell['branch_length'])
# plt.show()


# # estimate Ne (selPop Ne should be lower due to selection; maximise with only outlier trees)
# estimate.maximize(glike_fun, glike_x0, bounds=glike_bounds)


### HARMONIC SEARCH ESTIMATION
# from datetime import datetime


# nTimeVars=4

# def change(vals):
#   x1 = {"t1":vals[0],"tadmix":vals[1],"tsplit":vals[2], "tend":vals[3], "NeA":vals[4], "NeB":vals[5], "NeC":vals[6], "admixture":vals[6]}
#   return x1

# def normalize(val, lim):
#     ret = []
#     for i in range(len(val)):
#       nval = 0
#       if lim[i][1] != lim[i][0]:
#         nval = (val[i] - lim[i][0]) / (lim[i][1] - lim[i][0])
#       ret.append(nval)
#     return ret


# def denormalize(val, lim):
#     ret = []
#     for i in range(len(val)):
#         nval = (val[i]) * (lim[i][1] - lim[i][0]) + lim[i][0]
#         if (nval < lim[i][0]):
#             nval = lim[i][0]
#         elif (nval > lim[i][1]):
#             nval = lim[i][1]
#         ret.append(nval)

#     ret[0] = int(ret[0])
#     ret[1] = int(ret[1])
#     ret[2] = int(ret[2])
#     return ret

# def prob_function(pos):
#   vals = denormalize(pos[:], limits[:])
#   xt = change(vals)
#   #print("checking for",str(xt))
#   x, logp = estimate.estimateSingle(glike_fun,xt,bounds=glike_bounds)
#   return logp

# limits = glike_bounds[:]



# def harmonic(pHMCR,pPAR,pBW,filename,maxMin):
#   maxTime = maxMin*60
#   f = open(filename, "w")
#   print("simulating ihs", filename,"with",maxTime, flush=True)
#   starttime = datetime.now()
#   fname = "orighar"+str(rNum)+".txt"
  
#   listHarmonics = hsearch.loadhms(fname, limits,nTimeVars)[:]
  
#   for i in range(len(listHarmonics)):
#     # print("calculating for harmonic",i)
#     listHarmonics[i].score = prob_function(listHarmonics[i].values)
#     #print("calculating original harmoic",i,flush=True)

#   print(str(len(listHarmonics))+" original harmonics done",flush=True)
#   listHarmonics = hsearch.msorted(listHarmonics)[:]
#   """
#     j = 0
#      for i in listHarmonics:
#          print("For harnomic",j,"=",i.values,"score",i.score)
#          j += 1
#     """

#   bniter, his = hsearch.oloadhistory(fname)
#   # print("Before iterations ",bniter)
#   # print("prior history",his)

#   #rNum = int(input("Enter run folder"))
#   newharfile = "nprocess/newhar" + str(expId) + "exprun" + str(rNum) + "par"+ str(pPAR)+".txt"
#   hsearch.updateFile(newharfile, bniter, his, listHarmonics)
#   cbest = listHarmonics[-1]
#   cdiff = (datetime.now() - starttime).total_seconds()
#   print("#HAR#", str(cbest.values), flush=True)
#   print("#HAR#", cbest.score, flush=True)
#   print("#HAR#", cdiff / 60, flush=True)
#   f.write(str(cdiff / 60) + "#" + str(cbest.score) + "#" + str(cbest.values) + "\n")
#   print("cdiff initial",cdiff,flush=True)
#   ppTime = cdiff
  
#   if cdiff >= maxTime:
#      return

#   cutTimes = []
#   i = cdiff + 60
#   while i <= maxTime:
#     cutTimes.append(i)
#     i += 60
    
#   i = 0
#   j = 0
#   print(cutTimes,"cut times",flush=True)
#   optIter = (cutTimes[-1]-ppTime)/(ppTime/70)
#   print("predicted iterations",optIter,flush=True)
#   nsiter = 500000000
 


#   while i < nsiter:
#     fname = newharfile
    
#     listHarmonics = hsearch.loadhms(fname, limits,nTimeVars)[:]
#     listHarmonics = hsearch.msorted(listHarmonics)[:]
#     bniter, his = hsearch.loadhistory(fname)

#     vals = hsearch.newHarmonic(listHarmonics, limits, pHMCR, hsearch.ihsPar(pPAR, bniter,optIter), hsearch.ihsBw(pBW, bniter,optIter),nTimeVars)
#     # vals = hsearch.newHarmonic(listHarmonics, limits, 0.9, 0.3, 0.01)
#     ## print(vals)
#     # print(hsearch.denormalize(vals,limits))
#     prob = prob_function(vals)
#     # print("has probability",prob)
#     listHarmonics = hsearch.updatehms(listHarmonics, vals, prob)[:]
#     his.append(str(listHarmonics[-1].score) + ",0")
#     hsearch.updateFile(newharfile, bniter + 1, his, listHarmonics)
#     # print("updated file",i,"iteration")
#     i += 1

#     cbest = listHarmonics[-1]
#     cdiff = (datetime.now() - starttime).total_seconds()
#     #print(cdiff,"iteration",i,flush=True)
#     if cdiff > cutTimes[j]:
#       print("#HAR#", str(cbest.values), flush=True)
#       print("#HAR#", cbest.score, flush=True)
#       print("#HAR#", cdiff / 60, flush=True)
#       # f.write(str(cdiff / 60) + "," + str(cbest.score) + "\n")
#       # f.write(str(cdiff / 60) + "," + str(cbest.values) + "\n")
#       f.write(str(cdiff / 60) + "#" + str(cbest.score) + "#" + str(cbest.values) + "\n")
#       #writeMore(listHarmonics,fileident,cdiff)
      
#       j += 1
      
#       if j >= len(cutTimes):
#         print("iterations done =", i)
#         f.close()
#         print("finished ihs",flush=True)
#         return

#   f.close()
#   print("should not reach here?",flush=True)
#   print("finished ihs",flush=True)
#   rNum=1
#   expId=1
#   inpKappa=10000
#   runMinh=40
  
# hsearch.createOrigHarmonics("orighar"+str(1)+".txt",8,[6],1,75)
# harmonic(0.9,0.35,0.1,"crun"+str(1)+"/"+"exp"+str(expId)+"k"+str(inpKappa)+"har35.txt",runMinh)    

    



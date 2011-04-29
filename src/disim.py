# -*- coding: utf-8 -*-
#!/usr/bin/env python
r"""
Diffusion of innovation simulation.

This simulation models the processes discussed by Abrahamson and Rosenkopf in
[AR1997]_ and [RA1999]_.

:Author: Christopher Kirkos
:Date: 04/19/2011
:Version: |release|


Setup Model Procedure
---------------------

Choose turbulence (te) and complexity (ce) for environment
Choose distribution parameters based on te and ce

TODO List
---------

TODO: Incorporate reputational variance. Fn of Nx position or random selection
(1999 model [RA1999]_)
Bi,k = Ii + (A * sum(ri * Di,k-1)) for A >= 0
ri = reputation of org i
Di,k = 1 if org i has adopted by cycle k; 0 otherwise

TODO: Have global (A) and fixed (Ai) ambiguity experiments ([RA1999]_ p. 367)

TODO: Incorporate profitability assessment (1999 model [RA1999]_)

Implementation Details
----------------------
"""

from __future__ import division
from graphgen import generateARCorePeriph, drawAdoptionNetworkGV
from random import shuffle, choice, gauss
from itertools import product, cycle
from os.path import exists, join as pathjoin
from os import mkdir
import csv
from data import Data
from pylab import plt, rcParams

# Set global matplotlib style parameters
rcParams['legend.fontsize'] = 10
rcParams['lines.markersize'] = 3

# 1997 model: 3 sets of simulations:
#  1. Basic model of faddish diffusion
#  2. Relax assumption of equal sensitivity to information creating bandwagon
#     pressure
#  3. Model based on learning instead of fads

def run1997ThresholdModel(trickleDirection="down"):
	"""Runs the initial threshold model	from [AR1997]_"""
	# number of nodes in network
	numberOfNodes = 31
	# ratio of the number of nodes in the core to nodes in the periphery
	cpRatio = 1/3
	numCoreNodes = int(round(numberOfNodes*cpRatio))
	
	# set base location for output
	outFilePath = "/home/prima/Development/tmp/disim/out7"
	if not exists(outFilePath):
		mkdir(outFilePath)
	
	# set `drawNetworkImages` to true to output _ALL_ resulting DOT and PNG files
	drawNetworkImages = False
	
	# Record the results of every trial as a record in a CSV file 
	# Fields/columns (ordered): 
	#   # periphery ties, Ai, trial #, # core adopters, total # core nodes, 
	#   # periph adopters, total # periph nodes
	expTrialLogOutfile = "experimentTrialLog-n%d.csv" % numberOfNodes
	expTrialLogFileP = file(pathjoin(outFilePath,expTrialLogOutfile), "w")
	expTrialLogCSV = csv.writer(expTrialLogFileP)

	# Keep track of the average diffusion and density of the multiple trials
	# for each level of ambiguity Ai. Keep this in memory to graph later.
	# {ai : (avg peripheral diffusion, avg peripheral density), ... }
	experimentCaseLog = {}
		
	# Save trial data to csv too, stream to file
	# Columns: Ai, avg peripheral diffusion, avg peripheral density,
	# avg core diffusion 
	expCaseLogOutfile = "experimentCaseLog-n%d.csv" % numberOfNodes
	expCaseLogOutfileP = file(pathjoin(outFilePath,expCaseLogOutfile), "w")
	expCaseLogCSV = csv.writer(expCaseLogOutfileP)
	
	# total undirected simple edges = (n*(n-1))/2
	# total peripheral ties = ties between periphery, 
	# and ties between periph and core
	# == total graph ties - total ties between core only 
	# totalPossiblePeriphTies = (len(periphNodes)*(len(periphNodes)-1))/2
	totalPossibleTies = (numberOfNodes*(numberOfNodes-1))/2
	totalPossibleCoreTies = (numCoreNodes*(numCoreNodes-1))/2
	totalPossiblePeriphTies = totalPossibleTies - totalPossibleCoreTies
					
	# "We permitted the number of these ties to vary from 0 to 185 in intervals
	# of 5." ([AR1997]_ pp. 297-298)
	# IE. peripheryTies_i = xrange(0,185, 5) # Original use by authors
	# Instead, we scale this with the number of peripheral nodes so we can vary
	# the Network size if we want to.
	peripheryTies_i = xrange(0,int(totalPossiblePeriphTies),5)
	# TODO: parameterize the interval, currently set static to '5'
	
	# "In this first simulation, A_i was fixed to the same value for all
	# firms, but this value was permitted to vary between 1 and 5 in
	# intervals of 1." ([AR1997]_ p. 298)	
	A_i = xrange(1,6)	# [1, 2, 3, 4, 5]
	
	# generate all combinations of the # of ties and Ai for experimentation
	cases = product(peripheryTies_i, A_i)
	
	# "For each case, we ran 100 trials and calculated the average number of 
	# adopters in the focal and non-focal strata" ([AR1997]_ p. 298)
	trials = 100
	
	# A case is a combination of the number
	for pties,Ai in cases:
		if not experimentCaseLog.has_key(Ai):
			# default lists for experiment log.
			# peripheral diffusion, peripheral density, core diffusion
			experimentCaseLog[Ai] = ([],[],[])
		
		peripheralDiffusion = Data()
		peripheralDensity = Data()
		coreDiffusion = Data()
		
		# Generate a new network for each case
		Gorig = generateARCorePeriph(numCoreNodes,
									numberOfNodes-numCoreNodes, pties)
		trial = 1
		while trial<=trials:
			# make a copy of the generated graph b/c the simulation modifies 
			# the graph
			G = Gorig.copy()
			# "Assessed profits were drawn randomly from a normal distribution
			# with mean -1.0 and standard deviation 1.0" ([AR1997]_ p. 298)
			mu,sigma = -1.0,1.0
			
			# set the assessed profit for each node from normal distribution, I_i,
			# and the weight of bandwagon pressure, A_i.
			for a in G.nodes():
				G.node[a]['I'] = gauss(mu,sigma)
				G.node[a]['A'] = Ai
				G.node[a]['adopted'] = False
				G.node[a]['influence'] = []
			
			coreNodes = [a for a in G.nodes() if G.node[a]['core']]
			periphNodes = [a for a in G.nodes() if not G.node[a]['core']]
			
			if trickleDirection == "down":			
				# select a random core node for trickle-down diffusion
				seedNode = choice(coreNodes)
			else:
				seedNode = choice(periphNodes)
			
			G.node[seedNode]['adopted'] = True
			
			# start simulation
			while True:
				# we only need to evaluate agents that have not yet adopted
				agents = [n for n in G.nodes() if not G.node[n]['adopted']]
				
				# TODO: Option for simultaneous updating vs incremental.
				# Make a temp copy of the graph here, so that we have a
				# snapshot of the last round. Then, when complete, replace the 
				# simulation graph with the newly updated graph.
				# This would simulate simultaneous updating, instead of  
				# incremental.
				
				# Uniform random agent activation
				# TODO: Find out how activation occurred in the AR1997 model.
				shuffle(agents)
				madeChange = False
				for a in agents:
					# will agent a adopt?
					# compute B_i,k = I_i + (A_i * P_k-1)
					neighbors = G.neighbors(a)
					adoptedNeighbors = [n for n in neighbors \
										if G.node[n]['adopted'] == True]
					# In the 1997 fad model, Pk1 is the number of neighbor  
					# adopters divided by the total number of agents in the  
					# network (potential adopters)
					Pk1 = len(adoptedNeighbors)/G.number_of_nodes()
					Bik = G.node[a]['I'] + (Ai * Pk1)
					if Bik > 0:
						# Adopt if Bik was assessed > 0
						G.node[a]['adopted'] = True
						G.node[a]['influence'] = adoptedNeighbors
						madeChange = True
				
				# stop after no more agents can be influenced by bandwagon
				if not madeChange: 
					break
			
			if drawNetworkImages:
				# save resulting graph image to file
				outImgFilename = "n%d-PTies%d-Ai%d-Trial%d" % (numberOfNodes, 
															pties, Ai, trial)
				drawAdoptionNetworkGV(G, 
						writeFile=pathjoin(outFilePath, outImgFilename+".dot"),
						writePng=pathjoin(outFilePath, outImgFilename+".png"))
			
			# compute adopters in focal and non-focal strata
			numCoreAdopters = len([a for a in coreNodes 
									if G.node[a]['adopted']])
			periphNodes = [a for a in G.nodes() if not G.node[a]['core']]
			numPeriphAdopters = len([a for a in periphNodes  
									if G.node[a]['adopted']])
			# record experiment results
			expTrialLogCSV.writerow([pties, Ai, trial, numCoreAdopters, 
									len(coreNodes), numPeriphAdopters, 
									len(periphNodes)])
			
			peripheralDiffusion.addDatum(numPeriphAdopters/len(periphNodes))
			peripheralDensity.addDatum(pties/totalPossiblePeriphTies)
			coreDiffusion.addDatum(numCoreAdopters/len(coreNodes))
			trial += 1
		
		# TODO: Save experiment case log to file instead of in-memory only		
		experimentCaseLog[Ai][0].append(peripheralDensity.average)
		experimentCaseLog[Ai][1].append(peripheralDiffusion.average)
		experimentCaseLog[Ai][2].append(coreDiffusion.average)
		
		expCaseLogCSV.writerow((Ai, peripheralDensity.average,
									peripheralDiffusion.average,
									coreDiffusion.average))
	
	expTrialLogFileP.close()
	expCaseLogOutfileP.close()
	
	# Generate plot of Peripheral Diffusion vs. Nx density beyond the core
	# x axis: pties/total possible ties
	# y axis: # peripheral adopters / # periphery nodes
	
	mc = marker_cycle()
	fig = plt.figure()
	ax = fig.add_subplot(111)
	for Ai in experimentCaseLog.keys():
		# x axis is peripheral density, y axis is the peripheral diffusion
		x,y = experimentCaseLog[Ai][0], experimentCaseLog[Ai][1]
		ax.plot(x,y, label="Ambiguity=%d"%Ai, marker=mc.next())
	
	ax.set_xlabel("Network Density Beyond the Core")
	ax.set_ylabel("Peripheral Diffusion")
	ax.legend(loc="best")
	ax.set_title("Extent of Peripheral Diffusion for Varying Ambiguity and "
				"Network Density\n(Averaged over %d trials)" % trials)
	
	outPlotFilename = "DensityPlot-n%d.png" % (numberOfNodes)
	fig.savefig(pathjoin(outFilePath, outPlotFilename))
	

def run1997NxModel():
	"""Runs the modified network model from [AR1997]_"""
	pass

def run1999Model():
	"""Runs the expanded network model from [RA1999]_"""
	pass

def marker_cycle():
	""" Return an infinite, cycling iterator over the available marker 
	symbols.
	
	This is wrapped in a function to make sure that you get a new iterator
	that starts at the beginning every time you request one. This function is
	meant for use with Matplotlib.
	"""
	return cycle([
        'o','^','s','D','p','d','+','x','1','2','3','4','h',
        'H','|','_','v','<','>'])

if __name__ == "__main__":
	run1997ThresholdModel()


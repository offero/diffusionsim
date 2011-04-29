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
from plotting import createCoreDiffusionPlot, createPeripheralDiffusionPlot

from random import shuffle, choice, gauss
from itertools import product
from os.path import exists, join as pathjoin
from os import makedirs
import csv
from data import Data
from collections import defaultdict

# 1997 model: 3 sets of simulations:
#  1. Basic model of faddish diffusion
#  2. Relax assumption of equal sensitivity to information creating bandwagon
#     pressure
#  3. Model based on learning instead of fads

def run1997ThresholdModel(trickleDirection="down", numberOfNodes=31,
						outFilePath="/home/prima/Development/tmp/disim/out"):
	"""Runs the initial threshold model	from [AR1997]_"""
	# ratio of the number of nodes in the core to nodes in the periphery
	cpRatio = 1/3
	numCoreNodes = int(round(numberOfNodes*cpRatio))
	
	if not exists(outFilePath):
		makedirs(outFilePath)
	
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
	# {ai : (avg peripheral diffusion, avg peripheral density, 
	#        avg core diffusion), ... }
	experimentCaseLog = defaultdict(lambda: [[],[],[]])
		
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
	# IE. peripheryTies_i = xrange(0,185, 5) # For original use by the authors
	# Instead, we scale this with the number of peripheral nodes so we can vary
	# the Network size if we want to. We will also get the entire density range,
	# which is more computationaly expensive.
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
	
	expCaseLogOutfileP.close()
	expTrialLogFileP.close()
	
	periphDiffPlotTitle = "Extent of Peripheral Diffusion for Varying Ambiguity"\
				" and Network Density\n(Averaged over %d trials)" % trials
	createPeripheralDiffusionPlot(experimentCaseLog, outFilePath, 
								periphDiffPlotTitle)
	
	coreDiffPlotTitle = "Extent of Core Diffusion for Varying Ambiguity"\
				" and Network Density\n(Averaged over %d trials)" % trials
	createCoreDiffusionPlot(experimentCaseLog, outFilePath, 
								coreDiffPlotTitle)
	

def loadCaseLog(expCaseLogOutfilePath):
	"""Regenerate experiment case log structure from output log file."""
	
	# load data from output file
	expCaseLogOutfileP = file(expCaseLogOutfilePath, "r+")
	expCaseLogCSV = csv.reader(expCaseLogOutfileP)
	# Columns: Ai, avg peripheral diffusion, avg peripheral density,
	# avg core diffusion 
	
	experimentCaseLog = defaultdict(lambda: [[],[],[]])
	
	for Ai, pdiff, pdens, cdiff in expCaseLogCSV:		
		experimentCaseLog[Ai][0].append(pdiff)
		experimentCaseLog[Ai][1].append(pdens)
		experimentCaseLog[Ai][2].append(cdiff)
	
	return experimentCaseLog



from optparse import OptionParser, make_option
from sys import argv

def parseCommandLine():
	"""
	Available commands:
		simulate 
			-d, --direction=up/down/both
			-n, --nodes=<integer>
			-dots, --output-dot-files
			-pngs, --output-png-files
		plotstats 
			-i, --input-file=caseLogFile.csv
		plotnetwork 
			-i, --input-file=dotfile.dot
	
	Global options:
		Output directory: -o --output-dir
	
	"""
	
	command = argv[1]
	assert(command in ("simulate", "plotstats", "plotnetwork"))
	
	optlist=[
		# for simulate:
    	make_option("-d", "--direction", type="choice", 
					choices=("up","down","both"), dest="direction",
					default="down",
                 	help="Diffusion direction, either 'up' or 'down'."),
		make_option("-n", "--nodes", type="int", dest="numberOfNodes", 
					default=31,
					help="Output to directory. Default is the current working"\
						 " directory."),
		make_option("-dots", "--output-dot-files", action="store_true",
					dest="dots", default=False, 
					help="Output networks as .dot files for Graphviz"),
		make_option("-pngs", "--output-png-files", action="store_true",
					dest="pngs", default=False, 
					help="Generate Graphviz visualization of networks."),
		# for plotstats and plotnetwork
		make_option("-i", "--input-file", type="string", dest="inputFile", 
					help="Input file."),
		# for all
    	make_option("-o", "--output-dir", type="string", dest="outputDir", 
					default=".",
					help="Output to directory. Default is the current working"\
						 " directory."),		
	]
	parser = OptionParser(option_list=optlist)
	
	(options, args) = parser.parse_args()
	
	if command == "simulate":		
		trickleDirections = [options.direction,]
		if options.direction == "both":
			trickleDirections=["up","down"]
			
		for td in trickleDirections:
			run1997ThresholdModel(trickleDirection=td, 
					outFilePath=options.outputDir+"Trickle-%s-Simulation"%td)
	
	if command == "plotstats":
		experimentCaseLog = loadCaseLog(options.inputFile)
		createPeripheralDiffusionPlot(experimentCaseLog, options.outputDir)
		createCoreDiffusionPlot(experimentCaseLog, options.outputDir)
	
	if command == "plotnetwork":
		# TODO: Implement network plot
		pass


if __name__ == "__main__":
	parseCommandLine()



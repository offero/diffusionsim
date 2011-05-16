# -*- coding: utf-8 -*-
#!/usr/bin/env python
#
# Copyright (C) 2011 Christopher Kirkos. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

r"""
Diffusion of innovation simulation.

This simulation models the processes discussed by Abrahamson and Rosenkopf in
[AR1997]_ and [RA1999]_.

:Author: Christopher Kirkos
:Date: 04/19/2011
:Version: |release|

Implementation
--------------
"""

from __future__ import division

from graphgen import generateARCorePeriph, drawAdoptionNetworkGV
from plotting import createCoreDiffusionPlot, createPeripheralDiffusionPlot
from stats import possibleTies, runOLSRegression1997
from graphsearch import findWeaknessesAndPressurePoints, GRAPH_FILTERS,\
						clearWPPCache

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
						trials=100, cpRatio=1/3,
						outFilePath="/home/prima/Development/tmp/disim/out",
						dots="none", pngs="none"):
	"""Runs the initial threshold model	from [AR1997]_
	
	:param str trickleDirection: The direction of trickle simulation. This
									decides whether the seed adopter is in the
									core or periphery for the simulation.
	:param int numberOfNodes: The number of nodes in the generated network.
	:param int trials: The number of trials to run for each unique set of 
						   parameters (Periphery ties, Ai).
	:param float cpRatio: Ratio of the number of nodes in the core to nodes 
						  in the periphery.
	:param str outFilePath: The path where output files and sub-directories 
							   are created.
	:param str dots: Condition to output DOT files of the influence networks.
					 Possible values are "all", "wpp", and "none". "all" 
					 outputs files for each trial, "wpp" for only trials
					 whose resulting graphs have boundary conditions, and 
					 "none" for no outut.
	:param str pngs: Condition to output PNG files of the influence networks
					 (same values/conditions as *dots* argument).
	
	.. note::
		"For each case, we ran 100 trials and calculated the average number of 
		adopters in the focal and non-focal strata" ([AR1997]_ p. 298)		
	"""
	# Determine number of core nodes 
	numCoreNodes = int(round(numberOfNodes*cpRatio))
	
	targetSegment = 'periphery' if trickleDirection=="down" else "core"
	
	if not exists(outFilePath):
		makedirs(outFilePath)
	
	dotFilter = GRAPH_FILTERS[dots](targetSegment=targetSegment)
	pngFilter = GRAPH_FILTERS[pngs](targetSegment=targetSegment)
	
	# ***** The Experiment Trial Log *****
	# Record the results of every trial as a record in a CSV file 
	# Fields/columns (ordered): 
	# (0) # periphery ties, (1) Ai, (2) trial #, (3) # core adopters, 
	# (4) total # core nodes, (5) # periph adopters, (6) total # periph nodes,
	# (7) # boundary weaknesses, (8) # boundary pressure points
	expTrialLogOutfile = "experimentTrialLog-n%d.csv" % numberOfNodes
	expTrialLogFileP = file(pathjoin(outFilePath,expTrialLogOutfile), "w")
	expTrialLogCSV = csv.writer(expTrialLogFileP)
	# ************************************

	# ***** The Experiment Case Log *****
	# Keeps track of the average diffusion and density of the multiple trials
	# for each level of ambiguity Ai. Keep this in memory to graph later.
	# {ai : (avg peripheral diffusion, avg peripheral density, 
	#        avg core diffusion), ... }
	experimentCaseLog = defaultdict(lambda: [[],[],[]])		
	# Save trial data to csv too, stream to file
	# Columns: (0) Ai, (1) avg peripheral diffusion, (2) avg peripheral density,
	# (3) avg core diffusion 
	expCaseLogOutfile = "experimentCaseLog-n%d.csv" % numberOfNodes
	expCaseLogOutfileP = file(pathjoin(outFilePath,expCaseLogOutfile), "w")
	expCaseLogCSV = csv.writer(expCaseLogOutfileP)	
	# Note: The Experiment Case Log is primarily used to generate the  
	# peripheral/core diffusion graphs.	
	# ************************************
	
	# For calculating peripheral density
	totalPossiblePeriphTies = possibleTies(numberOfNodes, numCoreNodes)[2]
					
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
			
			# set the assessed profit (I_i) for each node from normal 
			# distribution, and the weight of bandwagon pressure (A_i).
			for a in G.nodes():
				G.node[a]['I'] = gauss(mu,sigma)
				G.node[a]['A'] = Ai
				G.node[a]['adopted'] = False
				G.node[a]['influence'] = []
			
			coreNodes = [a for a in G.nodes() if 'core' in \
										G.node[a]['segments']]
			periphNodes = [a for a in G.nodes() if 'core' not in \
												G.node[a]['segments']]
			
			# select a random core node as an adopter for trickle-down 
			# diffusion or a random peripheral node for trickle-up diffusion
			seedNode = choice(coreNodes) if trickleDirection == "down" \
										 else choice(periphNodes)
			G.node[seedNode]['adopted'] = True
			
			# Start simulation
			while True:
				# Only evaluate agents that have not yet adopted
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
					Bik = G.node[a]['I'] + (G.node[a]['A'] * Pk1)
					if Bik > 0:
						# Adopt if Bik was assessed > 0
						G.node[a]['adopted'] = True
						G.node[a]['influence'] = adoptedNeighbors
						madeChange = True
				
				# stop after no more agents can be influenced by bandwagon
				if not madeChange: 
					break
			
			# Find the boundary weaknesses and pressure points
			weaknesses, ppoints = findWeaknessesAndPressurePoints(G, 
						targetSegment="periphery" if trickleDirection=="down" \
						else "core")
			
			if pngs != "none" or dots != "none":
				# save resulting graph image to file
				outImgFilename = "n%d-PTies%d-Ai%d-Trial%d" % (numberOfNodes, 
															pties, Ai, trial)
				writeFileDot = pathjoin(outFilePath, outImgFilename+".dot") \
								if dotFilter(G) else None
				writeFilePng = pathjoin(outFilePath, outImgFilename+".png") \
								if pngFilter(G) else None
				drawAdoptionNetworkGV(G, 
									  writeFile=writeFileDot,
									  writePng=writeFilePng)
			
			# compute adopters in focal and non-focal strata
			numCoreAdopters = len([a for a in coreNodes 
									if G.node[a]['adopted']])
			numPeriphAdopters = len([a for a in periphNodes  
									if G.node[a]['adopted']])
			# record experiment results
			expTrialLogCSV.writerow([pties, Ai, trial, numCoreAdopters, 
									len(coreNodes), numPeriphAdopters, 
									len(periphNodes), len(weaknesses), 
									len(ppoints)])
			
			peripheralDiffusion.addDatum(numPeriphAdopters/len(periphNodes))
			peripheralDensity.addDatum(pties/totalPossiblePeriphTies)
			coreDiffusion.addDatum(numCoreAdopters/len(coreNodes))
			
			clearWPPCache() # resources about this graph are no longer needed
			trial += 1
		
		experimentCaseLog[Ai][0].append(peripheralDensity.average)
		experimentCaseLog[Ai][1].append(peripheralDiffusion.average)
		experimentCaseLog[Ai][2].append(coreDiffusion.average)
		
		expCaseLogCSV.writerow((Ai, peripheralDensity.average,
									peripheralDiffusion.average,
									coreDiffusion.average))
	
	expCaseLogOutfileP.close()
	expTrialLogFileP.close()
	
	periphDiffPlotTitle = "Extent of Peripheral Diffusion for Varying "\
				"Ambiguity and Network Density\n(Averaged over %d trials)"\
				 % trials
	createPeripheralDiffusionPlot(experimentCaseLog, outFilePath, 
								periphDiffPlotTitle)
	
	coreDiffPlotTitle = "Extent of Core Diffusion for Varying Ambiguity"\
				" and Network Density\n(Averaged over %d trials)" % trials
	createCoreDiffusionPlot(experimentCaseLog, outFilePath, 
								coreDiffPlotTitle)
	
	fullRegressionAnalysis(outFilePath, expTrialLogOutfile, trickleDirection)
	


def fullRegressionAnalysis(outFilePath, expTrialLogOutfile, trickleDirection):
	"""Perform regression analysis on all the identified combinations of values
	from the [AR1997]_ paper (With/without boundary conditions, all/high/low
	network densities, full/limited to 185 peripheral ties).
	
	This function should effectively run the necessary regression analyses 
	to regenerate all analysis tables presented in the paper.
	"""
	
	pTieRanges = (None, (0,185))
	densityRanges = (None, (0,0.5), (0.5,1))
	boundaryConds = (True,False)
	conditionCombos = product(pTieRanges, densityRanges, boundaryConds)
	
	for cond in conditionCombos:
		runOLSRegression1997(pathjoin(outFilePath, expTrialLogOutfile),
						trickleDirection=trickleDirection,
						peripheralTieRange = cond[0],
						densityRange = cond[1],
						withBoundaryAnalysis=cond[2],
                        outFilePath=outFilePath)
		

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
#from sys import argv

def parseCommandLine():
	"""
	Available commands:
		simulate 
			-d, --direction=up/down/both
			-n, --nodes=<integer>
			-t, --trials=<integer>
			-D, --dots=all,wpp
			-P, --pngs=all,wpp
		plotstats 
			-i, --input-file=caseLogFile.csv
		plotnetwork 
			-i, --input-file=dotfile.dot
	
	Global options:
		Output directory: -o --output-dir
	
	
	`simulate` runs the simulation
	`plotstats` takes a case log file (CSV) and produces a graph file (PNG)
	`plotnetwork` takes a DOT file and produces a network visualization (PNG)
	"""
	
	graphOutputChoices = GRAPH_FILTERS.keys() # from graphsearch module
	
	optlist=[
		# for simulate commands:
    	make_option("-d", "--direction", type="choice", 
					choices=("up","down","both"), dest="direction",
					default="down",
                 	help="Diffusion direction, either 'up' or 'down'."),
		make_option("-n", "--nodes", type="int", dest="numberOfNodes", 
					default=31,
					help="Number of nodes in network."),
		make_option("-t", "--trials", type="int", dest="trials", 
					default=100,
					help="Number of times to run each network simulation per "\
					     "set of parameters"),
		make_option("-D", "--dots", 
					dest="dots", type="choice", choices=graphOutputChoices,
					default="none", 
					help="Output networks as .dot files for Graphviz" \
					"Possible values are 'all' and 'wpp'. 'wpp' filters "\
					"so that only graphs with weaknesses and pressure points "\
					"are generated." ),
					#action="store_true",
		make_option("-P", "--pngs",
					dest="pngs", type="choice", choices=graphOutputChoices,
					default="none",
					help="Generate Graphviz visualization of networks." \
					"Possible values are 'all' and 'wpp'."),
		# for plotstats and plotnetwork commands:
		make_option("-i", "--input-file", type="string", dest="inputFile", 
					help="Input file."),
		# for all commands:
    	make_option("-o", "--output-dir", type="string", dest="outputDir", 
					default=".",
					help="Base output directory. Default is the current "\
						 "working directory."),
	]
	usage = "usage: %prog <command> [options] \n\n"\
			"Command is one of 'simulate' or 'plotstats'."
	parser = OptionParser(option_list=optlist, usage=usage)
	
	(options, args) = parser.parse_args()
	
	assert(len(args)>0 and args[0] in ("simulate", "plotstats"))
	assert(options.dots in graphOutputChoices and \
			options.pngs in graphOutputChoices)
	
	command = args[0]
	
	if command == "simulate":		
		trickleDirections = [options.direction,]
		if options.direction == "both":
			trickleDirections=["up","down"]
		
		for td in trickleDirections:
			outputFilePath = pathjoin(options.outputDir,
									"Trickle-%s-Simulation"%td)
			run1997ThresholdModel(trickleDirection=td, 
					numberOfNodes=options.numberOfNodes,
					trials = options.trials, 
					outFilePath=outputFilePath,
					dots=options.dots, pngs=options.pngs)
	
	if command == "plotstats":
		experimentCaseLog = loadCaseLog(options.inputFile)
		createPeripheralDiffusionPlot(experimentCaseLog, options.outputDir)
		createCoreDiffusionPlot(experimentCaseLog, options.outputDir)
	
	if command == "plotnetwork":
		# TODO: Implement network plot. If not useful, remove feature.
		pass

	if command == "regress":
		# TODO: Implement regression command line options
		pass


if __name__ == "__main__":
	parseCommandLine()



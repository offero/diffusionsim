# -*- coding: utf-8 -*-
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

'''
Created on Apr 29, 2011

Plotting for DISim.

:Author: Christopher Kirkos
'''

from __future__ import division

from pylab import plt, rcParams
from os.path import join as pathjoin
from itertools import cycle

# Set global matplotlib style parameters
rcParams['legend.fontsize'] = 10
rcParams['lines.markersize'] = 3

def createPeripheralDiffusionPlot(experimentCaseLog, outFilePath,
                                 plotTitle=None):
    """Function to generate the Peripheral Diffusion vs Peripheral Density
    plot.
    """
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
    if plotTitle == None:
        ax.set_title("Extent of Peripheral Diffusion for Varying Ambiguity "
                     "and Network Density")
    else:
        ax.set_title(plotTitle)
    
    outPlotFilename = "Plot-PeripheralDiffusionVsDensity.png"
    fig.savefig(pathjoin(outFilePath, outPlotFilename))

def createCoreDiffusionPlot(experimentCaseLog, outFilePath, plotTitle=None):
    """Function to generate the Core Diffusion vs Peripheral Density
    plot. (Basically a clone of `createPeripheralDiffusionPlot`)
    """
    # Generate plot of Core Diffusion vs. Nx density beyond the core
    # x axis: pties/total possible ties
    # y axis: # core adopters / # core nodes
    
    mc = marker_cycle()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for Ai in experimentCaseLog.keys():
        # x axis is peripheral density, y axis is the core diffusion
        x,y = experimentCaseLog[Ai][0], experimentCaseLog[Ai][2]
        ax.plot(x,y, label="Ambiguity=%d"%Ai, marker=mc.next())
    
    ax.set_xlabel("Network Density Beyond the Core")
    ax.set_ylabel("Core Diffusion")
    ax.legend(loc="best")
    if plotTitle == None:
        ax.set_title("Extent of Core Diffusion for Varying Ambiguity "
                     "and Network Density")
    else:
        ax.set_title(plotTitle)
    
    outPlotFilename = "Plot-CoreDiffusionVsDensity.png"
    fig.savefig(pathjoin(outFilePath, outPlotFilename))
    

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
    
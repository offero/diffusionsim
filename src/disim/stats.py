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
Created on May 2, 2011

Statistic generation for DISim

:Author: Christopher Kirkos
'''

from __future__ import division

import numpy as np
import scikits.statsmodels as sm
import numexpr as ne
#from itertools import product
from os.path import isdir, join as pathjoin

def standardizeCoeff(A, sample=True):
    """Element-wise, subtract the mean of the values from each value and 
    divide by the standard deviation.
    
    :param numpy.Array A: A numpy array.
    :param bool sample: Whether the `Array` `A` represents a sample or 
                           population (True if sample).
    :returns: A copy/view of the input Array with values standardized.
    :rtype: numpy.Array
    """
    sp = 1 if sample else 0 # sample standard deviation or population?
    return (A-np.mean(A))/np.std(A,ddof=sp)

def optimizedMinmax(A, newMin=0.0, newMax=1.0):
    """An optimized version of `minmax` using `numexpr`."""
    mn, mx = np.min(A), np.max(A)
    return ne.evaluate("(((A-mn)/(mx-mn))*(newMax-newMin))+newMin")

def minmax(A, newMin=0.0, newMax=1.0):
    """Min-Max normalization.
    
    :param numpy.Array A: A 1D numpy Array object.
    :param float newMin: The new minimum value.
    :param float newMax: The new maximum value.
    :returns: A copy/view of the input Array with values normalized.
    :rtype: numpy.Array
    """
    mn, mx = np.min(A), np.max(A)
    return (((A-mn)/(mx-mn))*(newMax-newMin))+newMin

def runOLSRegression1997(expTrialLogFilePath,
                         peripheralTieRange = (0,185),
                         densityRange = None,
                         withBoundaryAnalysis=False,
                         outFilePath=None):
    """Perform Ordinary Least Squares (OLS) regression on the trial output
    data to determine which parameters have the most effect on the extent
    of peripheral diffusion.
    
    :param str expTrialLogFilePath: The full path to the experiment's trial
                                       log file.
    :param boolean withBoundaryAnalysis: Whether or not to include the boundary
                                         weaknesses and pressure points in the
                                         multiple regression.
    :param tuple peripheralTieRange: Limit the regression analysis to the 
                                     records where the number of peripheral
                                     ties is within (including) the specified
                                     range (min, max). If `None` is specified, 
                                     then the range defaults to the entire 
                                     range of the sample (no restriction).
    :param tuple densityRange: Limit regression analysis to the records where
                               the peripheral density is within (including)
                               the specified range (Eg. *(0,0.5)*). If `None` 
                               is specified, then no density restriction is 
                               placed on the records. 
    
    .. note:: 
        The original authors only simulated with the number of peripheral
        ties within range [0,185]. In their second set of simulations, they 
        performed regressions on the same simulation but split the data set 
        by density of peripheral ties. They split it in half with the first 
        set having density <=0.5 and the second set with records having 
        density >0.5. This is the reasoning behind the `peripheralTieRange` 
        and `densityRange` function parameters.
    
    """
    # 7 Columns (from log file):
    # (0) periphery ties 
    # (1) ambiguity Ai
    # (2) trial #
    # (3) # core adopters
    # (4) total # core nodes 
    # (5) # periph adopters
    # (6) total # periph nodes
    # (7) # boundary weaknesses
    # (8) # boundary pressure points
    
    # Test output path, create new output name
    if outFilePath and isdir(outFilePath):
        ptr = peripheralTieRange
        ptiesStr = "{0:d}To{1:d}PTies-".format(ptr[0],ptr[1]) if ptr else ""
        dr = densityRange 
        pdensStr = "PDensity{0:.2f}To{1:.2f}-".format(dr[0],dr[1]) if dr else ""
        boundStr = "{0:s}Boundaries".format( "With" if withBoundaryAnalysis \
                                         else "Without" )
        fname = "Regression-{0:s}{1:s}{2:s}.txt".format(ptiesStr,pdensStr,boundStr)
        outFilePath = pathjoin(outFilePath,fname)
        
    # Import trial data from CSV file
    trialLogArray = np.genfromtxt(expTrialLogFilePath, dtype=np.float, 
                                  delimiter=',')
    
    # y is the dependent variable, the number of peripheral adopters
    y = yorig = trialLogArray[:,5]
    
    x1 = trialLogArray[:,1] # Ambiguity
    #x1 = optimizedMinmax(x1, newMin=-1.0, newMax=0.0) # normalize ambiguity
    
    x2 = trialLogArray[:,3]/trialLogArray[:,4] # Core diffusion
    # Core diffusion = core adopters / core nodes
    #   -> element-wise division on arrays
    
    # numexpr optimized version to calculate peripheral density
    x3 = optimizedCalcPeriphDensity(trialLogArray[:,0], trialLogArray[:,4], 
                                    trialLogArray[:,6])
    
    yname = "Peripheral diffusion"
    xnames = ["Ambiguity", "Core diff.", "Per. dens."]
    
    indep = [x1,x2,x3] # independent variables    
    if withBoundaryAnalysis:
        x4 = trialLogArray[:,7] # weaknesses
        x5 = trialLogArray[:,8] # pressure points
        indep.extend([x4,x5])
        xnames.extend(["Weaknesses", "Press. Pnts"])
    
    # Begin down selecting records based on input conditions
    pties = trialLogArray[:,0]
    ptr0,ptr1 = peripheralTieRange if peripheralTieRange else \
                                                (np.min(pties), np.max(pties)) 
    # selecting on pties, the # peripheral ties
    tieMask = np.ma.masked_inside(pties, ptr0, ptr1)
    
    densA = optimizedCalcNxDensity(trialLogArray[:,0], trialLogArray[:,4], 
                           trialLogArray[:,6])
    dr0,dr1 = densityRange if densityRange else (np.min(densA), np.max(densA))
    # Todo: ask authors which density they used? Nx or peripheral..
    # (originally) selecting on x3, the peripheral density independent variable
    # selecting on densA, the network density
    densityMask = np.ma.masked_inside(densA, dr0, dr1)
    
    combinedMaskVals = tieMask.mask & densityMask.mask

    indep = map(lambda x: x[combinedMaskVals], indep)
    y = y[combinedMaskVals]
    yorig = yorig[combinedMaskVals]
    
    # ensure we didn't down select to 0 sized arrays
    for x in indep + [y,yorig]:
        if len(x) <= 0:
            return None
    
    # standardize all coefficients
    indep = map(standardizeCoeff, indep)
    y = standardizeCoeff(y)
    
    # X represents the independent control variables
    X = np.array(indep, dtype=np.float).transpose()
    
    olsRegression = sm.OLS(y,X)
    olsFit = olsRegression.fit()
    regstats = (olsFit, (np.mean(yorig), np.std(yorig), np.min(yorig), 
                         np.max(yorig)))
    
    if outFilePath != None:
        with file(outFilePath, 'w') as outFileP:
            outFileP.write("Regression Summary\n")
            outFileP.write(olsFit.summary(yname=yname, xname=xnames))
            outFileP.write("\nMean: %f\nStdDev: %f\nMin: %f\nMax: %f\n" % \
                           (regstats[1]) )
    
    return regstats

def possibleTies(numberOfNodes, numCoreNodes):
    """ Calculates the max number of ties for a network with given node 
    statistics.
    
    Total undirected simple edges = (n*(n-1))/2
    Total peripheral ties = ties between nodes in the periphery and themselves,
    AND ties between nodes in the periphery and in the core core.
    
    Total peripheral ties is calculated as the ... 
    Total possible ties in the network minus the total ties between core nodes.
    
    :param integer numberOfNodes: The number of nodes in the network.
    :param integer numCoreNodes: The number of nodes in the core.
    :return: Tuple of the number of total possible ties in
             the network, the total possible core ties,
             and the total possible peripheral ties.
             
    .. todo:: 
        Write latex set notation to represent this calculation in the 
        docs.
    """
    totalPossibleTies = int( (numberOfNodes*(numberOfNodes-1))/2 )
    totalPossibleCoreTies = int( (numCoreNodes*(numCoreNodes-1))/2 )
    totalPossiblePeriphTies = int( totalPossibleTies - totalPossibleCoreTies )
    return (totalPossibleTies, totalPossibleCoreTies, totalPossiblePeriphTies)

def calcPerpipheralDensity(pties,coreNodes,periphNodes):
    """Calculate peripheral density, which is the number of ties in the
    periphery divided by the number of possible ties in the periphery. 
    
    All input values are scalars, not arrays/matrices.
    
    :param integer pties: The number of ties in the periphery.
    :param integer coreNodes: The number of nodes in the core.
    :param integer periphNodes: The number of nodes in the periphery.
    """
    totalNodes = coreNodes+periphNodes
    totalPossiblePeriphTies = possibleTies(totalNodes, coreNodes)[2]        
    return pties/totalPossiblePeriphTies
    
def optimizedCalcPeriphDensity(pt,cn,pn):
    """An optimized version of `calcPeripheralDensity` using the numexpr
    module.
    
    All input parameters are numpy Arrays.
    """
    tn = cn+pn # total nodes
    # totalPossibleTies = (tn*(tn-1))/2.0
    # totalPossibleCoreTies = (cn*(cn-1))/2.0
    # tot possible periph ties = totalPossibleTies - totalPossibleCoreTies
    # actual peripheral ties / tot possible periph ties  
    # pt/(( (tn*(tn-1))/2.0 ) - ( (cn*(cn-1))/2.0 ))
    return ne.evaluate("pt/( ((tn*(tn-1))/2.0) - ((cn*(cn-1))/2.0) )")


def optimizedCalcNxDensity(pt,cn,pn):
    "An optimized version of `calcNxDensity`"
    tn = cn+pn
    return ne.evaluate("(pt+((cn*(cn-1))/2.0)) / ((tn*(tn-1))/2.0)")

def calcNxDensity(pties,coreNodes,periphNodes):
    """Calculate the density of a network with the given number of peripheral
    ties, core nodes, and peripheral nodes"""
    totalPossibleTies, totalPossibleCoreTies, totalPossiblePeriphTies = \
                                possibleTies(coreNodes+periphNodes, coreNodes)
    return (pties+totalPossibleCoreTies)/totalPossibleTies

if __name__ == "__main__":
    import disim
    basepath = '/home/prima/Development/tmp/disim/testreg'
    infile = 'experimentTrialLog-n31.csv'
    disim.fullRegressionAnalysis(basepath, infile)
    

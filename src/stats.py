'''
Created on May 2, 2011

Statistic generation for DISim

:Author: Christopher Kirkos
'''

from __future__ import division

import numpy as np
import scikits.statsmodels as sm
import numexpr as ne

def standardizeCoeff(A, sample=True):
    """Element-wise, subtract the mean of the values from each value and 
    divide by the standard deviation.
    
    :param numpy.Array A: A numpy array.
    :param bool sample: Whether the `Array` `A` represents a sample or 
                           population (True if sample).
    :return numpy.Array: A copy/view of the input Array with values 
                         standardized.
    """
    sp = 1 if sample else 0 # sample standard deviation or population?
    return (A-np.mean(A))/np.std(A,ddof=sp)

def optimizedMinmax(A, newMin=0.0, newMax=1.0):
    """An optimized version of `minmax`_ using `numexpr`."""
    mn, mx = np.min(A), np.max(A)
    return ne.evaluate("(((A-mn)/(mx-mn))*(newMax-newMin))+newMin")

def minmax(A, newMin=0.0, newMax=1.0):
    """Min-Max normalization.
    
    :param numpy.Array A: A 1D numpy Array object.
    :param float newMin: The new minimum value.
    :param float newMax: The new maximum value.
    :return numpy.Array: A copy/view of the input Array with values normalized.
    """
    mn, mx = np.min(A), np.max(A)
    return (((A-mn)/(mx-mn))*(newMax-newMin))+newMin

def runOLSRegression1997(expTrialLogFilePath, tieLimit=None, 
                         withBoundaryAnalysis=False):
    """Perform Ordinary Least Squares (OLS) regression on the trial output
    data to determine which parameters have the most effect on the extent
    of peripheral diffusion.
    
    :param str expTrialLogFilePath: The full path to the experiment's trial
                                       log file.
    :param boolean withBoundaryAnalysis: Whether or not to include the boundary
                                         weaknesses and pressure points in the
                                         multiple regression.
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
    
    # Import trial data from CSV file
    trialLogArray = np.genfromtxt(expTrialLogFilePath, dtype=np.float, 
                                  delimiter=',')
    
    if tieLimit != None:
        trialLogArray = trialLogArray[:tieLimit]
    
    # y is the dependent variable, the number of peripheral adopters
    y = trialLogArray[:,5]
    
    # y is the dependent variable, peripheral diffusion (adopters/nodes)
    #y = trialLogArray[:,5]/trialLogArray[:,6]
    
    x1 = trialLogArray[:,1] # Ambiguity
    #x1 = optimizedMinmax(x1, newMin=-1.0, newMax=0.0) # normalize ambiguity
    
    x2 = trialLogArray[:,3]/trialLogArray[:,4] # Core diffusion
    # Core diffusion = core adopters / core nodes
    #   -> element-wise division on arrays
    
    # numexpr optimized version to calculate peripheral density
    x3 = optimizedCalcPeriphDensity(trialLogArray[:,0], trialLogArray[:,4], 
                                    trialLogArray[:,6])
    
    # numpy only version of calculating peripheral density
    #vCalcPDensity = np.vectorize(calcPerpipheralDensity)    
    #x3 = vCalcPDensity(trialLogArray[:,0], trialLogArray[:,4], 
    #                   trialLogArray[:,6]) # Peripheral density
    
    y,x1,x2,x3 = [standardizeCoeff(a) for a in [y,x1,x2,x3]]
    
    indep = [x1,x2,x3] # independent variables    
    if withBoundaryAnalysis:
        x4 = standardizeCoeff(trialLogArray[:,7]) # weaknesses
        x5 = standardizeCoeff(trialLogArray[:,8]) # pressure points
        indep.extend([x4,x5])
    
    # X represents the independent control variables
    X = np.array(indep, dtype=np.float).transpose()
    
    olsRegression = sm.OLS(y,X)
    olsFit = olsRegression.fit()
    yp = trialLogArray[:,5] # the original dependent variable, y
    stats = (olsFit, (np.mean(yp), np.std(yp), np.min(yp), np.max(yp)))
    return stats

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
    :return Tuple(int,int,int): Tuple of the number of total possible ties in
                                the network, the total possible core ties,
                                and the total possible peripheral ties.
    """
    # TODO: Write latex set notation to represent this calculation in the docs.
    totalPossibleTies = (numberOfNodes*(numberOfNodes-1))/2
    totalPossibleCoreTies = (numCoreNodes*(numCoreNodes-1))/2
    totalPossiblePeriphTies = totalPossibleTies - totalPossibleCoreTies
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
    """An optimized version of `calcPeripheralDensity`_ using the numexpr
    module.
    
    All input parameters are numpy Arrays.
    """
    tn = cn+pn # total nodes
    # totalPossibleTies = (tn*(tn-1))/2.0
    # totalPossibleCoreTies = (cn*(cn-1))/2.0
    # tot possible periph ties = totalPossibleTies - totalPossibleCoreTies
    # actual peripheral ties / tot possible periph ties  
    # pt/(( (tn*(tn-1))/2.0 ) - ( (cn*(cn-1))/2.0 ))
    return ne.evaluate("pt/(( (tn*(tn-1))/2.0 ) - ( (cn*(cn-1))/2.0 ))")

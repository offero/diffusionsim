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
'''
Created on May 9, 2011

:author: Christopher Kirkos

About
-----

This module tests the graph generation capabilities of DISim.

Testing `disim.graphgen.DICorePeriphNxGenerator`_
+++++++++++++++++++++++++++++++++++++++++++++++++ 

DICorePeriphNxGenerator is an infinite generator, to test it we need to supply
it with all combinations of boundary values for each parameter.

Important values are negative, zero, positive but < max nodes/edges, 
max nodes/edges, and > max nodes/edges. 

Test implementation
+++++++++++++++++++
'''
from nose.tools import raises
from itertools import permutations, product
from disim.graphgen import DICorePeriphNxGenerator
from disim.stats import possibleTies

@raises(Exception)
def createInvalidGraph(numCore, numPeriph, pties):
    nxGen = DICorePeriphNxGenerator(numCore, numPeriph, pties)


def createAndTestValidGraph(numCore, numPeriph, pties):
    nxGen = DICorePeriphNxGenerator(numCore, numPeriph, pties)
    totNxTies,totCoreTies,totPeriphTies = possibleTies(numCore+numPeriph, 
                                                       numCore)
    # test the generator for at least 2 networks
    for i in range(2):
        G = nxGen.next()
        # Check for total number of nodes  
        assert(G.number_of_nodes() == numCore+numPeriph)
        # Check the number of edges created
        assert(G.number_of_edges() == totCoreTies+pties)
        assert(G.number_of_edges() <= totNxTies)
        
        periphNodes = [n for n in G.nodes() if 'core' not in \
                                            G.node[n]['segments']]        
        # Check for number of peripheral edges beyond core created
        assert(len(G.edges(periphNodes)) == pties)
        assert(len(G.edges(periphNodes)) <= totPeriphTies)
        # Check number of peripheral nodes
        assert(len(periphNodes) == numPeriph)
        
        coreNodes = [n for n in G.nodes() if 'core' in G.node[n]['segments']]
        # Check number of core nodes
        assert(len(coreNodes) == numCore)
        

def testGraphGeneration():
    numCoreNodesL = numPeriphNodesL = [1, 2, 10, 100]
    ptiesL = [0,1]
    errorCombosNeg1 = permutations([-1,10,10], 3) # -1 in each parameter
    
    #okCombosZero = permutations([0,10,10], 3) # 0 in each parameter
    #okCombosOne = permutations([1,10,10], 3) # 1 in each parameter 
    
    for numCore,numPeriph,pties in errorCombosNeg1:
        yield createInvalidGraph, numCore, numPeriph, pties
    
    # Asking for edges without nodes
    yield createInvalidGraph, 0, 0, 1
    yield createInvalidGraph, 0, 0, 10
    # Case: When 0 for either core or periph and >0 for peripheral ties
    #       This case is invalid, as no peripheral ties are possible 
    yield createInvalidGraph, 1, 0, 1
    yield createInvalidGraph, 0, 1, 1
    yield createInvalidGraph, 1, 0, 1
    yield createInvalidGraph, 0, 1, 10
    yield createInvalidGraph, 1, 0, 10
    
    # Though it is valid to ask for 0 peripheral ties...
    yield createAndTestValidGraph, 0, 0, 0
    yield createAndTestValidGraph, 1, 0, 0
    yield createAndTestValidGraph, 0, 1, 0
    
    for numCore,numPeriph in product(numCoreNodesL, numPeriphNodesL):
        #totalPossibleTies, totalPossibleCoreTies, totalPossiblePeriphTies        
        totNx,totCore,totPeriph = possibleTies(numCore+numPeriph, numCore)
        # test lower boundary of peripheral ties
        for pties in ptiesL:
            yield createAndTestValidGraph, numCore, numPeriph, pties
        # test upper boundary of peripheral ties
        yield createAndTestValidGraph, numCore, numPeriph, totPeriph
        # Trying to create more ties than possible should be an error...
        yield createInvalidGraph, numCore, numPeriph, totPeriph+1
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    
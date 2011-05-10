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

The purpose of this module is to test the networks search functions that
detect weaknesses and pressure points.
'''

from disim.graphsearch import findWeaknessesAndPressurePoints

import networkx as nx
from itertools import combinations

def testDetectWeaknessesAndPressurePoints():
    
    numCore = 4
    connectedCore = combinations(range(numCore), 2) # nodes 0,1,2,3 fully connected
    edges = [(0,5),(1,5),(2,5),(3,5),(5,6), # pressure point on 5
             (0,7), (2,7), (6,7), # weak and pressure point node 7
             (1,8), (6,8)             
             ]
    edges.extend(connectedCore)
    G = nx.Graph(data=edges)
    
    for node in G.nodes():
        if node < numCore:
            G.node[node]['segments'] = ['core',]
        else:
            G.node[node]['segments'] = ['periphery',]
        G.node[node]['I'] = -2.0
        G.node[node]['A'] = 3.0
    
    G.node[7]['I'] = 1.0 # weak node
    
    weaknesses,ppoints = findWeaknessesAndPressurePoints(G)
    #print weaknesses
    #print ppoints
    assert(weaknesses == [7])
    assert(ppoints == [5,7])
    
if __name__ == "__main__":
    testDetectWeaknessesAndPressurePoints()
    
    
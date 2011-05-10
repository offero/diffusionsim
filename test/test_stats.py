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

There are standard and optimized versions of calculations in the stats 
module. These tests are made to ensure that they compute as expected and that
the optimized versions match the output of the standard versions.

'''
from __future__ import division

from numpy import vectorize, array

from disim.stats import ( standardizeCoeff, 
                          calcNxDensity, 
                          calcPerpipheralDensity, 
                          optimizedCalcNxDensity, 
                          optimizedCalcPeriphDensity,
                          possibleTies )


def testNetworkDensity():
    
    vCalcNxDensity = vectorize(calcNxDensity)
    vCalcPerpipheralDensity = vectorize(calcPerpipheralDensity)
        
    def dens(pt,cn,pn):
        tn = cn+pn
        return (((cn*(cn-1))/2)+pt) / ((tn*(tn-1))/2)
    
    # pties, coreNodes, periphNodes
    data = ((10,10,20),
            (10,10,30),
            (10,10,40),
            (1,1,1),
            (1,10,10),
            (10,1,1)
           )
    
    # manually calculate network density for each example in data
    c10p10dens = array([dens(pt,cn,pn) for pt,cn,pn in data], dtype=float)
    
    # create a numpy array that will be given to the function
    dataA = array(data, dtype=float).T
    res1 = optimizedCalcNxDensity(dataA[0],dataA[1],dataA[2])
    res2 = vCalcNxDensity(dataA[0],dataA[1],dataA[2])
    r1eq2 = res1 == res2
    r2eqc = res2 == c10p10dens
    r1eqc = res1 == c10p10dens
    assert (False not in r1eq2)
    assert (False not in r2eqc)
    assert (False not in r1eqc)
    
    # Compare the standard and optimized versions of peripheral density
    pdensS = vCalcPerpipheralDensity(dataA[0],dataA[1],dataA[2])
    pdensO = optimizedCalcPeriphDensity(dataA[0],dataA[1],dataA[2])
    pdensSvO = pdensS == pdensO
    assert (False not in pdensSvO)
    
    
if __name__ == "__main__":
    testNetworkDensity()         
    
    
    
    
    
    
    
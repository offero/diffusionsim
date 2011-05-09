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

"""
A clone of Robert Axtell's Data.java.
"""

from __future__ import division
import sys, math

class Data(object):
    """An object that maintains running mean, variance/standard deviation, min, 
    max, and number of samples.
    """
    def __init__(self):
        self.N = 0
        self.min = sys.float_info.max
        self.max = sys.float_info.min
        self.sum = 0.0
        self.sum2 = 0.0

    def addDatum(self, val):
        "Add a datum, updating the running statistics."
        self.N = self.N + 1
        if val < self.min: self.min = val
        if val > self.max: self.max = val
        self.sum = self.sum + val
        self.sum2 = self.sum2 + (val**2)

    @property
    def average(self):
        "The average of the data (accessed as a property)."
        return self.sum/self.N if self.N > 0 else 0.0

    @property
    def variance(self):
        "The variance of the data (accessed as a property)."
        if self.N > 1:
            avg = self.average
            arg = self.sum2 - (self.N * (avg**2))
            return arg / (self.N - 1)
        else:
            return 0.0

    @property
    def stdDev(self):
        "The standard deviation of the data (accessed as a property)."
        return math.sqrt(self.variance)


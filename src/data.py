"""
A clone of Robert Axtell's Data.java class from assignments 2 and 3.
"""

from __future__ import division
import sys, math

class Data(object):
    def __init__(self):
        self.N = 0
        self.min = sys.float_info.max
        self.max = sys.float_info.min
        self.sum = 0.0
        self.sum2 = 0.0

    def addDatum(self, val):
        self.N = self.N + 1
        if val < self.min: self.min = val
        if val > self.max: self.max = val
        self.sum = self.sum + val
        self.sum2 = self.sum2 + (val**2)

    @property
    def average(self):
        return self.sum/self.N if self.N > 0 else 0.0

    @property
    def variance(self):
        if self.N > 1:
            avg = self.average
            arg = self.sum2 - (self.N * (avg**2))
            return arg / (self.N - 1)
        else:
            return 0.0

    @property
    def stdDev(self):
        return math.sqrt(self.variance)


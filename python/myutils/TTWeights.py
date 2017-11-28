#!/usr/bin/env python
import ROOT
import math

class TTWeights(object):

    def __init__(self):
        self.branches = [{'name': 'TTW', 'formula': self.getTTW}]

    def getBranches(self):
        return self.branches

    def getTTW(self, tree):
        TTW = 1
        if tree.nGenTop == 2:
            sf_top1 = math.exp(0.0615 - 0.0005*tree.GenTop_pt[0])
            sf_top2 = math.exp(0.0615 - 0.0005*tree.GenTop_pt[1])
            TTW = math.sqrt(sf_top1*sf_top2)
        return TTW

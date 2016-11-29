# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: reactions.py
@time: 11/29/2016 2:36 PM
"""

import species, math
import scipy.constants as sc

class reactions(object):

    def __init__(self, IS, TS, FS, T):
        self._IS = IS
        self._TS = TS
        self._FS = FS
        self._T = T


    @property
    def IS(self):
        return self._IS

    @IS.setter
    def IS(self,s):
        self._IS = s

    @property
    def FS(self):
        return self._FS

    @FS.setter
    def FS(self,s):
        self._FS = s

    @property
    def TS(self):
        return self._TS

    @TS.setter
    def TS(self,s):
        self._TS = s

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self,s):
        self._T = s

    @property
    def forwardE(self):
        return self._TS.energy - self._IS.energy

    @property
    def reverseE(self):
        return self._TS.energy - self._FS.energy

    @property
    def forwardK(self):
        return self.kHTST(self.forwardE)

    @property
    def reverseK(self):
        return self.kHTST(self.reverseE)

    def kHTST(self, E):
        k = sc.physical_constants['Boltzmann constant in eV/K'][0]
        h = sc.physical_constants['Planck constant in eV s'][0]
        return (k * self._T / h) * math.e ** ( - E / (k * self._T) )

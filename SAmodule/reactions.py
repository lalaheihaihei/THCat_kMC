# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: reactions.py
@time: 11/29/2016 2:36 PM
"""

import math

import scipy.constants as sc

from SAmodule import species


class reactions(object):

    def __init__(self, IS, TS, FS, T):
        self._IS = IS
        self._TS = TS
        self._FS = FS
        self._T = T
        self.__k = sc.physical_constants['Boltzmann constant in eV/K'][0]
        self.__h = sc.physical_constants['Planck constant in eV s'][0]
        self.__c = sc.physical_constants['speed of light in vacuum'][0]

    @property
    def IS(self):
        return self._IS

    @IS.setter
    def IS(self,s):
        if not isinstance(self._IS, species):
            raise ValueError('Error: input in reactions must be species.')
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
        return self._TS.correctEnergy - self._IS.correctEnergy

    @property
    def reverseE(self):
        return self._TS.correctEnergy - self._FS.correctEnergy

    @property
    def forwardK(self):
        return self.kHTST(self.forwardE, self._IS.freq, self._TS.freq)

    @property
    def reverseK(self):
        return self.kHTST(self.reverseE, self._FS.freq, self._TS.freq)

    def kHTST(self, E, freq1, freq2):

        productRatio = self.productRatio(freq1, freq2)
        #print("productRatio: ", productRatio, 'E: ', E)
        # See 'Fundamental ConCepts in Heterogeneous Catalysis' equation 4.47
        return (self.__k * self._T / self.__h) * productRatio * math.e ** ( - E / (self.__k * self._T) )

    def productRatio(self, freq1, freq2):
        product1 = product2 = 1
        for i in freq1:
            product1 *= (1 - math.e ** ((-1.23981e-4 * i) / (self.__k * self._T)))
        for i in freq2:
            product2 *= (1 - math.e ** ((-1.23981e-4 * i) / (self.__k * self._T)))
        #print(product1, product2)
        return product1/product2
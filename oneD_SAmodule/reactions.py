# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: reactions.py
@time: 11/29/2016 2:36 PM
"""

import math

import scipy.constants as sc

from SAmodule import species


class Reactions(object):

    def __init__(self, energies, reactions, temperature, freq):
        self._energies = energies
        self._reactions = reactions
        self._T = temperature
        self._freq = freq
        self.__k = sc.physical_constants['Boltzmann constant in eV/K'][0]
        self.__h = sc.physical_constants['Planck constant in eV s'][0]
        self.__c = sc.physical_constants['speed of light in vacuum'][0]

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self,s):
        self._T = s

    def ZPE(self):
        zpe = []
        for i in self._freq:
            zpe.append(sum(list(map(lambda x: 1/2 * x * 1.23981e-4, i))))
        return zpe


    @property
    def forwardE(self):
        zpe = self.ZPE()
        return (self._energies[1] + zpe[1]) - (self._energies[0] + zpe[0])

    @property
    def reverseE(self):
        zpe = self.ZPE()
        return (self._energies[1] + zpe[1]) - (self._energies[2] + zpe[2])

    @property
    def forwardK(self):
        return self.kHTST(self.forwardE, self._freq[0], self._freq[1])

    @property
    def reverseK(self):
        return self.kHTST(self.reverseE, self._freq[2], self._freq[1])

    def kHTST(self, E, freq1, freq2):
        productRatio = self.productRatio(freq1, freq2)
        # print("productRatio: ", productRatio, 'E: ', E)
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
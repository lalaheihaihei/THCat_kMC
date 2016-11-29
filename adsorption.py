# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: adsorption.py
@time: 11/29/2016 4:10 PM
a class of the
"""

import species, math, re
import scipy.constants as sc
import numpy as np

class adsorption(object):

    def __init__(self, IS, FS, T, P, gasKind = 'CO'):
        self._IS = IS
        self._FS = FS
        self._T = T
        self._P = P
        self._gasKind = gasKind


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
    def T(self):
        return self._T

    @T.setter
    def T(self,s):
        self._T = s

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self,s):
        self._P = s

    @property
    def gasKind(self):
        return self._gasKind

    @gasKind.setter
    def gasKind(self,s):
        self._gasKind = s


    def gas_mu(self):
        filename = 'JANAF_%s.txt' % (self._gasKind)
        s0 = h0 = h = s = mu_eV_P = 0
        with open(filename, 'r') as f:
            for line in f.readlines():
                line = line.split()
                #print(line)
                if int(line[0]) == 0:
                    s0 = float(line[1])
                    h0 = float(line[2])
                if int(line[0]) == int(self._T):
                    s = float(line[1]) - s0
                    h = float(line[2]) - h0
                mu_KjperMol = h - self._T * s / 1000 # KJ/mol
                mu_eV = mu_KjperMol * 0.010364272    # eV
                k = sc.physical_constants['Boltzmann constant in eV/K'][0]
                mu_eV_P = mu_eV + k * self._T * math.log(self._P/1.)
        return mu_eV_P


    def deltaG(self):
        #T_S_CO = -0.703269353 # at 300 K and 10e-3 bar
        isG = self._IS.energy + self.gas_mu()
        fsG = self._FS.energy
        return fsG - isG

    def adsorbK(self):
        mass = 0
        for i in self._gasKind:
            if i == 'C':
                mass += 12.01
            elif i == 'O':
                mass += 16.00
            elif i == '2':
                mass += 16.00
            else:
                raise ValueError('Please set molecular mass by yourself')
        s = 0.5  # sticking coefficient, we assume S = 0.5 for all the species
        A = (0.3 * 10e-9) ** 2  # the area of the adsorption site
        m = (mass / 1000) / sc.physical_constants['Avogadro constant'][0] # mass
        k = sc.physical_constants['Boltzmann constant'][0]  # boltzmann constant
        return ( s * self._P * A ) / math.sqrt( 2 * math.pi * m * k * self._T)

    def desorbK(self):
        k = sc.physical_constants['Boltzmann constant in eV/K'][0]
        K_equilibrium = math.e ** ( - self.deltaG() / ( k * self._T) ) \
            # The equilibrium constant K = k(forward)/k(backward)
        desorbK = self.adsorbK() / K_equilibrium
        return desorbK


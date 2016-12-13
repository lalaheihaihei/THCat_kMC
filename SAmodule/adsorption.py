# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: adsorption.py
@time: 11/29/2016 4:10 PM
a class of the
"""

import math
import scipy.constants as sc


class Adsorption(object):

    def __init__(self, ini_state, fin_state, temperature, pressure, gas_kind='CO'):
        self._IS = ini_state
        self._FS = fin_state
        self._T = temperature
        self._P = pressure
        self._gas_kind = gas_kind

    @property
    def ini_state(self):
        return self._IS

    @ini_state.setter
    def ini_state(self, s):
        self._IS = s

    @property
    def fin_state(self):
        return self._FS

    @fin_state.setter
    def fin_state(self, s):
        self._FS = s

    @property
    def temperature(self):
        return self._T

    @temperature.setter
    def temperature(self, s):
        self._T = s

    @property
    def pressure(self):
        return self._P

    @pressure.setter
    def pressure(self, s):
        self._P = s

    @property
    def gas_kind(self):
        return self._gas_kind

    @gas_kind.setter
    def gas_kind(self, s):
        self._gas_kind = s

    def gas_mu(self):
        filename = '../Thermo/JANAF_%s.txt' % self._gas_kind
        s0 = h0 = h = s = mu_ev_p = 0
        with open(filename, 'r') as f:
            for line in f.readlines():
                line = line.split()
                # print(line)
                if int(line[0]) == 0:
                    s0 = float(line[1])
                    h0 = float(line[2])
                if int(line[0]) == int(self._T):
                    s = float(line[1]) - s0
                    h = float(line[2]) - h0
                mu_kjpermol = h - self._T * s / 1000  # KJ/mol
                mu_ev = mu_kjpermol * 0.010364272     # eV
                k = sc.physical_constants['Boltzmann constant in eV/K'][0]
                mu_ev_p = mu_ev + k * self._T * math.log(self._P/1.)
        return mu_ev_p

    def delta_g(self):
        # T_S_CO = -0.703269353 # at 300 K and 10e-3 bar
        is_g = self._IS.energy + self.gas_mu()
        fs_g = self._FS.energy
        return fs_g - is_g

    def adsorb_k(self):
        if self._gas_kind == 'CO':
            mass = 28.01
        elif self._gas_kind == 'CO2':
            mass = 44.01
        elif self._gas_kind == 'O2':
            mass = 32.00
        elif self._gas_kind == 'H2':
            mass = 2.02
        elif self._gas_kind == 'N2':
            mass = 28.01
        elif self._gas_kind == 'NH3':
            mass = 17.03
        else:
            raise ValueError('Error, cannot recongnize the gas molecule.')
        s = 0.5  # sticking coefficient, we assume S = 0.5 for all the species
        area = (0.3 * 10e-9) ** 2  # the area of the adsorption site
        m = (mass / 1000) / sc.physical_constants['Avogadro constant'][0]  # mass
        k = sc.physical_constants['Boltzmann constant'][0]  # boltzmann constant
        return (s * self._P * area) / math.sqrt(2 * math.pi * m * k * self._T)

    def desorb_k(self):
        k = sc.physical_constants['Boltzmann constant in eV/K'][0]
        k_equilibrium = math.e ** (- self.delta_g() / (k * self._T)) \
            # The equilibrium constant K = k(forward)/k(backward)
        desorb_k = self.adsorb_k() / k_equilibrium
        return desorb_k

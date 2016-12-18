# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: adsorption.py
@time: 12/17/2016 20:32 PM
a class of the
"""


import math
import scipy.constants as sc
from oneD_SAmodule import adsorption


class Adsorption(object):

    def __init__(self, energies, reactions, temperature, pressure, gas_kind='CO'):
        self._energies = energies
        self._reactions = reactions
        self._T = temperature
        self._P = pressure
        self._gas_kind = gas_kind

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
                while '' in line:
                    line.remove('')
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
        return self._energies[1] - (self._energies[0] + self.gas_mu())

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
        elif self._gas_kind == 'C2H2':
            mass = 26.038
        elif self._gas_kind == 'C2H4':
            mass = 28.054
        elif self._gas_kind == 'C2H6':
            mass = 30.07
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


def get_ads_rate(kmc, i, T, P):
    # print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ads')
    ads = adsorption.Adsorption(kmc.Energies[i], kmc.reactions[i], T, P[kmc.reactions[i][0][-1].strip()], kmc.reactions[i][0][-1].strip())
    # print(ads.adsorb_k(), ads.desorb_k())
    return ads.adsorb_k(), ads.desorb_k()


def get_des_rate(kmc, i, T, P):
    energylist_swap = [kmc.Energies[i][1], kmc.Energies[i][0]]
    reactions_list_swap = [kmc.reactions[i][1], kmc.reactions[i][0]]
    # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    # print(energylist_swap, reactions_list_swap, T, P[kmc.reactions[i][1][-1].strip()], kmc.reactions[i][1][-1].strip())
    ads = adsorption.Adsorption(energylist_swap, reactions_list_swap, T, P[kmc.reactions[i][1][-1].strip()], kmc.reactions[i][1][-1].strip())
    return ads.desorb_k(), ads.adsorb_k()


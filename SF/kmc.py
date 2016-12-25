# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py.py
@time: 12/24/2016 23:44 PM
"""

from oneD_SAmodule import adsorption
from oneD_SAmodule import reactions
from oneD_SAmodule import loop
from IO import Config


def add_reaction(kmc, T, P, rate_const = [], rate_const_dict = {}):
    """
    Calculate rate constants of all the elementary reaction step.

    :param kmc: instance kmc include all input information
    :param T: temperature
    :param P: pressure
    :return: rate_const: list include all k.
             rate_const_dict: dict like {'reaction step': k, '-0': 8.96e-09, '6': 3.69e+33, '-6': 10.31,}
    """
    for i in range(len(kmc.reactions)):
        # calculate forwards (k1) and reverses (k2) rate constants.
        if kmc.reactionsKind[i] == 'ads':
            k1, k2 = adsorption.get_ads_rate(kmc, i, T, P)
        elif kmc.reactionsKind[i] == 'react':
            k1, k2 = reactions.get_react_rate(kmc, i, T)
        elif kmc.reactionsKind[i] == 'des':
            k1, k2 = adsorption.get_des_rate(kmc, i, T, P)
        else:
            raise ValueError('Error: check the reactionKind')
        # store the rate constants in list rate_const and rate_const_dict
        for j in [k1, k2]:
            rate_const.append(j)
            rate_const_dict[str(i)], rate_const_dict['-' + str(i)] = k1, k2
    return rate_const, rate_const_dict


def main():
    kmc = Config.Parameters()

    # get reaction rate constants list
    rate_const, rate_const_dict = add_reaction\
        (kmc, kmc.Temperature, {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))})

if __name__ == "__main__":
    main()
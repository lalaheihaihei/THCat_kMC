# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py.py
@time: 12/17/2016 4:08 PM
"""
import math
import random
import sys

from oneD_SAmodule import adsorption
from oneD_SAmodule import reactions
# from SAmodule import species
from IO import Config
# from SAmodule import lattice


def test_reaction_number(kmc):
    """
    Test the input reaction number. Whether the n(ads+des+rea) = n(reactions)

    :param kmc: Instance of Parameters
    :return: None
    """
    if len(kmc.reactionsKind) != len(kmc.reactions):
        raise ValueError('Error: Please check the reactions number and reactionsKind number')
    return None


def add_list(kmc):
    return None


def add_reaction(kmc, T, P, rate_const = []):
    """

    :param kmc:
    :param T:
    :param P:
    :return:
    """
    print('####################################################################')
    print(kmc.Energies)
    print(kmc.reactions)
    print(kmc.reactionsKind)
    print(kmc.Freq)
    for i in range(len(kmc.reactions)):
        if kmc.reactionsKind[i] == 'ads':
            k1, k2 = adsorption.get_ads_rate(kmc, i, T, P)
        elif kmc.reactionsKind[i] == 'react':
            k1, k2 = reactions.get_react_rate(kmc, i, T)
        elif kmc.reactionsKind[i] == 'des':
            k1, k2 = adsorption.get_des_rate(kmc, i, T, P)
        else:
            raise ValueError('Error: check the reactionKind')
        for j in [k1, k2]:
            rate_const.append(j)
    return rate_const


def initialize_lattice(kmc):
    """
    ATTENTION: this section should be define by users.
    :param kmc:
    :return:
    """
    lat = []
    for i in range(kmc.latticeSize[0]):
        if i % 5. < 0.1: # if i/5 = 0
            lat.append('1')
        else:
            lat.append('0')
    return lat


def initialize_num_of_avail_sites(kmc, lat):
    n_avail = {}
    for i in range(len(kmc.reactions)):
        n_avail[str(i)] = 0
        n_avail[str(-i)] = 0
    for i in range(len(kmc.reactions)):
        for j in lat:
            if j == kmc.reactions[i][0][0].strip():
                n_avail[str(i)] += 1
            if j == kmc.reactions[i][-1][0].strip():
                n_avail[str(-i)] += 1
    print(n_avail)
    return n_avail


def do_kmc_loop(rate_const, lat, num_of_avail_sites):
    return None


def main():
    kmc = Config.Parameters()

    T = kmc.Temperature                                               # get temperature /K
    P = {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))}  # get pressure    /Bar

    # Test the reaction number.
    test_reaction_number(kmc)

    add_list(kmc)

    rate_const = add_reaction(kmc, T, P)
    print('the list of rate constant is ', rate_const)
    for i in range(len(kmc.reactionsKind)):
        print("For step %d:\tk(forward) = %.3e,\tk(backward) = %.3e," % ( i, rate_const[i*2], rate_const[i*2 + 1]))

    lat = initialize_lattice(kmc)
    print(lat)

    '''num_of_avail_sites should be like: [[0,2,3,6],[1,5],[4,7]...]'''
    num_of_avail_sites = initialize_num_of_avail_sites(kmc, lat)
    print(num_of_avail_sites)

    #time, productNum, coverageList = \
    #    do_kmc_loop(rate_const, lat, num_of_avail_sites)



if __name__ == "__main__":
    main()
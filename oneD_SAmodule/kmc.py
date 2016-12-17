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


def get_ads_rate(kmc, i, T, P):
    # print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ads')
    ads = adsorption.Adsorption(kmc.Energies[i], kmc.reactions[i], T, P[kmc.reactions[i][0][-1].strip()], kmc.reactions[i][0][-1].strip())
    # print(ads.adsorb_k(), ads.desorb_k())
    return ads.adsorb_k(), ads.desorb_k()


def get_react_rate(kmc, i, T):
    order_of_the_react = kmc.reactionsKind[:i].count('react')  # to judge which three freq need to be input.
    input_freq = []
    for j in [order_of_the_react*3 + 0, order_of_the_react*3 + 1, order_of_the_react*3 + 2]:
        input_freq.append(kmc.Freq[j])
    # print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    # print(kmc.Energies[i], kmc.reactions[i], T, input_freq)
    react = reactions.Reactions(kmc.Energies[i], kmc.reactions[i], T, input_freq)
    # print("For step %d:\tk(forwards) = %.3e,\tk(reverse) = %.3e" % (i, react.forwardK, react.reverseK))
    return react.forwardK, react.reverseK


def get_des_rate(kmc, i, T, P):
    energylist_swap = [kmc.Energies[i][1], kmc.Energies[i][0]]
    reactions_list_swap = [kmc.reactions[i][1], kmc.reactions[i][0]]
    # print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    # print(energylist_swap, reactions_list_swap, T, P[kmc.reactions[i][1][-1].strip()], kmc.reactions[i][1][-1].strip())
    ads = adsorption.Adsorption(energylist_swap, reactions_list_swap, T, P[kmc.reactions[i][1][-1].strip()], kmc.reactions[i][1][-1].strip())
    return ads.desorb_k(), ads.adsorb_k()


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
            k1, k2 = get_ads_rate(kmc, i, T, P)
        elif kmc.reactionsKind[i] == 'react':
            k1, k2 = get_react_rate(kmc, i, T)
        elif kmc.reactionsKind[i] == 'des':
            k1, k2 = get_des_rate(kmc, i, T, P)
        else:
            raise ValueError('Error: check the reactionKind')
        for j in [k1, k2]:
            rate_const.append(j)
    return rate_const


def initialize_lattice(kmc):
    return None


def initialize_num_of_avail_sites(kmc, lat):
    return None


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

    #lat = initialize_lattice(kmc)

    '''num_of_avail_sites should be like: [[0,2,3,6],[1,5],[4,7]...]'''
    #num_of_avail_sites = initialize_num_of_avail_sites(kmc, lat)

    #time, productNum, coverageList = \
    #    do_kmc_loop(rate_const, lat, num_of_avail_sites)



if __name__ == "__main__":
    main()
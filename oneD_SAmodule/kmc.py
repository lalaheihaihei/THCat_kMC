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
from oneD_SAmodule import loop
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


def add_reaction(kmc, T, P, rate_const = [], rate_const_dict = {}):
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
            rate_const_dict[str(i)], rate_const_dict['-' + str(i)] = k1, k2
    return rate_const, rate_const_dict


def initialize_lattice(kmc):
    """
    ATTENTION: this section should be define by users.
    :param kmc:
    :return:
    """
    lat = []
    for i in range(kmc.latticeSize[0]):
        if i % (int(kmc.latticeSize[0])//float(kmc.num_SA)) < 0.1: # int(kmc.latticeSize[0])//float(kmc.num_SA) is the interval of SAs;
            lat.append('1')
        else:
            lat.append('0')
    lat[0] = '-1'  # non-periodic
    lat[-1] = '-1'  # non-periodic
    return lat


def count_of_forwards_reactions(kmc, n_avail, lat, i, j):
    """
    one list element in n_avail["-0"] means one avail site for this reaction. the element is also an element that
    contains reaction and site information.

    "0" : [[i, j, j], [i, j, j-1]] means the reaction "0" happens on site [j, j] and [j, j-1], which means
    one site reaction and two sites reaction, respectively.
    And the i store the reaction information as order in kmc.reactions

    :param kmc:
    :param n_avail:
    :param lat:
    :param i:
    :param j:
    :return:
    """
    # for case: 7 --> 8
    if len(kmc.reactions[i][0]) == 1 and kmc.reactions[i][0][0].strip() in kmc.kinds:
        if lat[j] == kmc.reactions[i][0][0].strip():
            n_avail[str(i)].append([str(i), j, j])
    # for case 7 + CO --> 8
    elif len(kmc.reactions[i][0]) == 2 and kmc.reactions[i][0][-1].strip() not in kmc.kinds:
        if lat[j] == kmc.reactions[i][0][0].strip():
            n_avail[str(i)].append([str(i), j, j])
    # for case 7 + 8 --> 9
    elif len(kmc.reactions[i][0]) == 2 and kmc.reactions[i][0][-1].strip() in kmc.kinds:
        if lat[j] == kmc.reactions[i][0][0].strip() and lat[j - 1] == kmc.reactions[i][0][1].strip():
            n_avail[str(i)].append([str(i), j, j-1])
        if lat[j] == kmc.reactions[i][0][0].strip() and lat[j + 1] == kmc.reactions[i][0][1].strip():
            n_avail[str(i)].append([str(i), j, j+1])
    # for case 7 + 8 + CO --> 10
    elif len(kmc.reactions[i][0]) == 3 and kmc.reactions[i][0][1].strip() in kmc.kinds:
        if lat[j] == kmc.reactions[i][0][0].strip() and lat[j - 1] == kmc.reactions[i][0][1].strip():
            n_avail[str(i)].append([str(i), j, j-1])
        if lat[j] == kmc.reactions[i][0][0].strip() and lat[j + 1] == kmc.reactions[i][0][1].strip():
            n_avail[str(i)].append([str(i), j, j+1])
    return n_avail


def count_of_reverse_reactions(kmc, n_avail, lat, i, j):
    # for case: 7 --> 8
    if len(kmc.reactions[i][1]) == 1 and kmc.reactions[i][1][0].strip() in kmc.kinds:
        if lat[j] == kmc.reactions[i][1][0].strip():
            n_avail['-'+str(i)].append(["-"+str(i), j, j])
    # for case 7 + CO --> 8
    elif len(kmc.reactions[i][1]) == 2 and kmc.reactions[i][1][-1].strip() not in kmc.kinds:
        if lat[j] == kmc.reactions[i][1][0].strip():
            n_avail['-'+str(i)].append(["-"+str(i), j, j])
    # for case 7 + 8 --> 9
    elif len(kmc.reactions[i][1]) == 2 and kmc.reactions[i][1][-1].strip() in kmc.kinds:
        if lat[j] == kmc.reactions[i][1][0].strip() and lat[j - 1] == kmc.reactions[i][1][1].strip():
            n_avail['-'+str(i)].append(["-"+str(i), j, j-1])
        if lat[j] == kmc.reactions[i][1][0].strip() and lat[j + 1] == kmc.reactions[i][1][1].strip():
            n_avail['-'+str(i)].append(["-"+str(i), j, j+1])
    # for case 7 + 8 + CO --> 10
    elif len(kmc.reactions[i][1]) == 3 and kmc.reactions[i][1][1].strip() in kmc.kinds:
        if lat[j] == kmc.reactions[i][1][0].strip() and lat[j - 1] == kmc.reactions[i][1][1].strip():
            n_avail['-'+str(i)].append(["-"+str(i), j, j-1])
        if lat[j] == kmc.reactions[i][1][0].strip() and lat[j + 1] == kmc.reactions[i][1][1].strip():
            n_avail['-'+str(i)].append(["-"+str(i), j, j+1])
    return n_avail


def initialize_num_of_avail_sites(kmc, lat):
    n_avail = {}                        # dict to record which react happen on which site {'-0':[1,5,8,...],'0':[]}
    for i in range(len(kmc.reactions)):
        n_avail[str(i)] = []            # list to record site
        n_avail['-'+str(i)] = []        # list to record site
    for i in range(len(kmc.reactions)): # loop all reactions
        for j in range(1, len(lat)-1):  # loop all lattice ,cutoff first and last to avoid exceeding of list range
            n_avail = count_of_forwards_reactions(kmc, n_avail, lat, i, j)
            n_avail = count_of_reverse_reactions(kmc, n_avail, lat, i, j)
    return n_avail


def main():
    kmc = Config.Parameters()

    T = kmc.Temperature                                               # get temperature /K
    P = {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))}  # get pressure    /Bar

    # Test the reaction number.
    test_reaction_number(kmc)

    add_list(kmc)

    rate_const, rate_const_dict = add_reaction(kmc, T, P)
    print('the list of rate constant is ', rate_const)
    for i in range(len(kmc.reactionsKind)):
        print("For step %d:\tk(forward) = %.3e,\tk(backward) = %.3e," % ( i, rate_const[i*2], rate_const[i*2 + 1]))

    lat = initialize_lattice(kmc)

    '''num_of_avail_sites should be like: [[0,2,3,6],[1,5],[4,7]...]'''
    num_of_avail_sites = initialize_num_of_avail_sites(kmc, lat)

    kmc_loop = loop.Loop(kmc, rate_const_dict, lat, num_of_avail_sites)
    print(kmc_loop.do_kmc_loop())



if __name__ == "__main__":
    main()
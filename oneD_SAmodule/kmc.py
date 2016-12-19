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


def add_reaction(kmc, T, P, rate_const = [], rate_const_dict = {}):
    """
    Calculate rate constants of all the elementary reaction step.
    :param kmc: object kmc include all input information
    :param T: temperature
    :param P: pressure
    :return: rate_const: list include all k.
             rate_const_dict: dict like {'reaction step': k, '-0': 8.96e-09, '6': 3.69e+33, '-6': 10.31,}
    """
    print('####################################################################')
    print("all elementary reactions energy are\n", kmc.Energies)
    print("all elementary reactions in list are\n:", kmc.reactions)
    # print(kmc.Freq)
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


def initialize_lattice(kmc, count_cut_num_of_active = 0):
    """
    ATTENTION: this section should be define by users.
    :param kmc:
    :return: the number_SA in config file does not the final active site num,
     because the code take first and last element in lat to be "-1" or PBC.
     Thus count_cut_num_of_active is used to record the cut number of active site,
     which is useful for final coverage calculation.
    """
    lat = []
    for i in range(kmc.latticeSize[0]):
        if i % (int(kmc.latticeSize[0])//float(kmc.num_SA)) < 0.1: # int(kmc.latticeSize[0])//float(kmc.num_SA) is the interval of SAs;
            lat.append('1')
        else:
            lat.append('0')
    for i in [0, -1]:
        if lat[i] == "1":
            count_cut_num_of_active += 1
    if kmc.periodic == "0":
        lat[0] = '-1'  # non-periodic
        lat[-1] = '-1'  # non-periodic
    elif kmc.periodic == "1":
        lat[0] = lat[-2]
        lat[-1] = lat[1]
    else:
        raise ValueError("Error: please check input periodic condition.")
    return lat, count_cut_num_of_active


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
        n_avail[str(i)] = []            # creat list to record site
        n_avail['-'+str(i)] = []        # creat list to record site
    for i in range(len(kmc.reactions)): # loop all reactions
        for j in range(1, len(lat)-1):  # loop all lattice ,cutoff first and last to avoid exceeding of list range
            # for periodic condition. lat[0] is lat[-2], lat[-1] is lat[1].
            # for nonperiodic conditon. lat[0] = lat[-1] = "-1"
            n_avail = count_of_forwards_reactions(kmc, n_avail, lat, i, j)
            n_avail = count_of_reverse_reactions(kmc, n_avail, lat, i, j)
    return n_avail


def main():
    kmc = Config.Parameters()

    T = kmc.Temperature                                               # get temperature /K
    P = {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))}  # get pressure    /Bar

    # Test the reaction number.
    test_reaction_number(kmc)

    rate_const, rate_const_dict = add_reaction(kmc, T, P)
    for i in range(len(kmc.reactionsKind)):
        print("For step %d:\tk(forward) = %.3e,\tk(backward) = %.3e," % ( i, rate_const[i*2], rate_const[i*2 + 1]))

    lat, count_cut_num_of_active = initialize_lattice(kmc)

    # num_of_avail_sites should be like: [[0,2,3,6],[1,5],[4,7]...]
    num_of_avail_sites = initialize_num_of_avail_sites(kmc, lat)

    kmc_loop = loop.Loop(kmc, rate_const_dict, lat, num_of_avail_sites, count_cut_num_of_active)
    print(kmc_loop.do_kmc_loop())



if __name__ == "__main__":
    main()
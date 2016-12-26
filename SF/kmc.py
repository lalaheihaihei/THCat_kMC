# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py.py
@time: 12/24/2016 23:44 PM
"""

from SF import adsorption
from SF import reactions
from SF import er
from SF import count_reaction
from SF import loop
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
        if kmc.reactions[i][1] == 'ads':
            k1, k2 = adsorption.get_ads_rate(kmc, i, T, P)
        elif kmc.reactions[i][1] == 'react':
            k1, k2 = reactions.get_react_rate(kmc, i, T)
        elif kmc.reactions[i][1] == 'des':
            k1, k2 = adsorption.get_des_rate(kmc, i, T, P)
        elif kmc.reactions[i][1] == 'er':
            k1, k2 = er.get_rate(kmc, i, T, P)
        else:
            raise ValueError('Error: check the reactionKind')

        # store the rate constants in list rate_const and rate_const_dict
        for j in [k1, k2]:
            rate_const.append(j)
            rate_const_dict[str(i)], rate_const_dict['-' + str(i)] = k1, k2
    return rate_const, rate_const_dict


def initialize_lattice(kmc, count_cut_num_of_active = 0):
    """
    ATTENTION: this section should be define by users.

    :param kmc: instance kmc include all input information
    :return: the number_SA in config file does not the final active site num,
     because the code take first and last element in lat to be "-1" or PBC.
     Thus count_cut_num_of_active is used to record the cut number of active site,
     which is useful for final coverage calculation.
    """
    lat = []
    if kmc.surface_kind == "rutile-110":
        for i in range(kmc.latticeSize[0]):
            lat.append([])
            for j in range(kmc.latticeSize[1]):
                if i % 2 == 0:
                    lat[i].append('18')
                elif i % 2 == 1:
                    lat[i].append('24')
        if kmc.periodic == "0":
            lat.insert(0, [])
            lat.insert(kmc.latticeSize[0] + 1, [])
            for j in range(kmc.latticeSize[1]):
                lat[0].append('-1')
                lat[-1].append('-1')
            for i in range(kmc.latticeSize[0] + 2):
                lat[i].insert(0,"-1")
                lat[i].insert(kmc.latticeSize[1] + 1,"-1")
        if kmc.periodic == "1":
            lat.insert(0, [])
            lat.insert(kmc.latticeSize[0] + 1, [])
            for j in range(kmc.latticeSize[1]):
                lat[0].append(lat[-2][j])
                lat[-1].append(lat[1][j])
            for i in range(kmc.latticeSize[0] + 2):
                lat[i].insert(0,lat[i][-1])
                lat[i].insert(kmc.latticeSize[1] + 1,lat[i][1])
    lat[3][3] = "1"
    return lat

def main():
    kmc = Config.Parameters()

    # get reaction rate constants list
    rate_const, rate_const_dict = add_reaction\
        (kmc, kmc.Temperature, {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))})

    # output rate constants
    for i in range(len(kmc.reactions)):
        print("For step %d: \t k(forward) = %.3e,\tk(backward) = %.3e," % (i+1, rate_const[i * 2], rate_const[i * 2 + 1]))

    # initialize lattice species
    lat= initialize_lattice(kmc)
    for i in lat:
        print(i)

    # get num_of_avail_sites that should be like: [[0,2,3,6],[1,5],[4,7]...]
    num_of_avail_sites = count_reaction.initialize_num_of_avail_sites(kmc, lat)
    for i in range(len(kmc.reactions)):
        print(str(i), '\t\t', num_of_avail_sites[str(i)])
        print('-'+str(i), '\t\t', num_of_avail_sites['-'+str(i)])

    # kmc_loop is an instance of class Loop in oneD_SAmodule.loop
    kmc_loop = loop.Loop(kmc, rate_const_dict, lat, num_of_avail_sites)
    # do loop
    kmc_loop.do_kmc_loop()
    print("###########################   normal termination   #####################################")
if __name__ == "__main__":
    main()
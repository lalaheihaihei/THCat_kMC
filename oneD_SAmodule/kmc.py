# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py.py
@time: 12/17/2016 4:08 PM
"""
import math
import random
import sys

from SAmodule import adsorption
from SAmodule import reactions
from SAmodule import species
from IO import Config
from SAmodule import lattice


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


def get_adsrate(kmc, i):
    return None


def get_reactrate(kmc, i):
    return None


def get_desrate(kmc, i):
    return None


def add_reaction(kmc, rate_const = []):
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
        if kmc.reactionsKind == 'ads':
            k1, k2 = get_adsrate(kmc, i)
        elif kmc.reactionsKind == 'react':
            k1, k2 = get_reactrate(kmc, i)
        elif kmc.reactionsKind == 'des':
            k1, k2 = get_desrate(kmc, i)
        else:
            raise ValueError('Error: check the reactionKind')
        (rate_const.append(k1)).append(k2)
        for j in [k1, k2]:
            rate_const.append(j)
    print('the list of rate constant is ', rate_const)
    return rate_const


def initialize_lattice(kmc):
    return None


def initialize_num_of_avail_sites(kmc, lat):
    return None


def do_kmc_loop(rate_const, lat, num_of_avail_sites):
    return None


def main():
    kmc = Config.Parameters()

    # T = kmc.Temperature                                               # get temperature /K
    # P = {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))}  # get pressure    /Bar

    # Test the reaction number.
    test_reaction_number(kmc)

    add_list(kmc)

    rate_const = add_reaction(kmc)

    lat = initialize_lattice(kmc)

    '''num_of_avail_sites should be like: [[0,2,3,6],[1,5],[4,7]...]'''
    num_of_avail_sites = initialize_num_of_avail_sites(kmc, lat)

    time, productNum, coverageList = \
        do_kmc_loop(rate_const, lat, num_of_avail_sites)



if __name__ == "__main__":
    main()
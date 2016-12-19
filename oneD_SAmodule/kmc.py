# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py.py
@time: 12/17/2016 4:08 PM
"""

from oneD_SAmodule import adsorption
from oneD_SAmodule import reactions
from oneD_SAmodule import loop
from IO import Config


def test_input_parameters(kmc):
    """
    Test the input parameters.
    (1) Whether the n(ads+des+rea) = n(reactions)
    (2) Whether the lattice is one dimension
    (3, 4) Check the reactions' number
    (5) Whether the corresponding energies and reaction equation for "react" is 3
    (6) Check the frequencies number

    :param kmc: Instance of Parameters
    :return: if no error, return None
    """
    if len(kmc.gas) != len(kmc.GasPressure):
        raise ValueError("Error: input gas number is not equal to GasPressure number")
    if kmc.latticeSize[1] != 1:
        raise ValueError("Error: oneD_SAmodule request latticeSize to be 1")
    if len(kmc.reactionsKind) != len(kmc.reactions):
        raise ValueError('Error: Please check the reactions number and reactionsKind number')
    if len(kmc.reactionsKind) != len(kmc.Energies):
        raise ValueError('Error: Please check the reactions number and Energies number')
    count_react = 0
    for i in range(len(kmc.reactionsKind)):
        if kmc.reactionsKind[i] == "react":
            if len(kmc.Energies[i]) != 3 or len(kmc.reactions[i]) != 3:
                raise ValueError("Error: Please check the 'react'")
            count_react += 1
    if count_react * 3 != len(kmc.Freq):
        raise ValueError("Error: the Frequency number is not right")
    print("#####################   input finish   #########################\n")
    return None


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
    for i in range(kmc.latticeSize[0]):
        # int(kmc.latticeSize[0])//float(kmc.num_SA) is the interval of SAs;
        # below "if" controls initial species, which need to be revised according to surface structure.
        if i % (int(kmc.latticeSize[0])//float(kmc.num_SA)) < 0.1:
            lat.append('1')
        else:
            lat.append('0')
    # WARNING: count_cut_num_of_active is used to cutoff active site in lat[0] and lat[-1]
    for i in [0, -1]:
        if lat[i] == "1":
            count_cut_num_of_active += 1
    # For non-periodic system, let the lat[0] and lat[-1] to be "-1",
    # and be careful that there should be no "-1" in kmc.kinds, thus the simulation will be cut at the head and tail.
    if kmc.periodic == "0":
        lat[0] = '-1'
        lat[-1] = '-1'
    # For periodic system, let the head and tail of lattice list repeat the "last" and "next" list.
    # see also update_lat(self) in loop.py
    elif kmc.periodic == "1":
        lat[0] = lat[-2]
        lat[-1] = lat[1]
    else:
        raise ValueError("Error: please check input periodic condition.")
    return lat, count_cut_num_of_active


def count_of_forwards_reactions(kmc, n_avail, lat, i, j):
    """
    this def append ['-6', 3, 3] values for keys in n_avail.

    one list element in n_avail["-0"] means one avail site for this reaction. the element is also an element that
    contains reaction and site information.

    "0" : [[i, j, j], [i, j, j-1]] means the reaction "0" happens on site [j, j] and [j, j-1], which means
    one site reaction and two sites reaction, respectively.
    And the i store the reaction information as order in kmc.reactions

    :param kmc: instance kmc include all input information
    :param n_avail: {'-6': [['-6', 3, 3]], '5': [], '-2': []...}
    :param lat: ['0', '0', '0', '1', '0', '0']
    :param i: loop of all reactions
    :param j: loop of all surface site. In oneD module, only one site or two site reaction are allowed.
    :return: n_avail
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
    """
    see notes in def count_of_forwards_reactions
    """
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
    """
    def to initialize number of available sites for specific reactions.
    :param kmc: instance kmc include all input information
    :param lat: ['0', '0', '0', '1', '0', '0']
    :return: n_avail dict like {'-6': [['-6', 3, 3]], '5': [], '-2': []...}
    """
    n_avail = {}                        # dict to record which react happen on which site {'-0':[1,5,8,...],'0':[]}
    for i in range(len(kmc.reactions)):
        n_avail[str(i)] = []            # creat list to record site
        n_avail['-'+str(i)] = []        # creat list to record site
    for i in range(len(kmc.reactions)): # loop all reactions
        for j in range(1, len(lat)-1):  # loop all lattice ,cutoff first and last to avoid double count in Periodic case
            n_avail = count_of_forwards_reactions(kmc, n_avail, lat, i, j)
            n_avail = count_of_reverse_reactions(kmc, n_avail, lat, i, j)
    return n_avail


def main():
    kmc = Config.Parameters()

    # Test the reaction number.
    test_input_parameters(kmc)

    # get reaction rate constants list
    rate_const, rate_const_dict = add_reaction\
        (kmc, kmc.Temperature, {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))})

    # output rate constants
    for i in range(len(kmc.reactionsKind)):
        print("For step %d:\tk(forward) = %.3e,\tk(backward) = %.3e," % (i, rate_const[i*2], rate_const[i*2 + 1]))

    # initialize lattice species
    lat, count_cut_num_of_active = initialize_lattice(kmc)

    # get num_of_avail_sites that should be like: [[0,2,3,6],[1,5],[4,7]...]
    num_of_avail_sites = initialize_num_of_avail_sites(kmc, lat)

    # kmc_loop is an instance of class Loop in oneD_SAmodule.loop
    kmc_loop = loop.Loop(kmc, rate_const_dict, lat, num_of_avail_sites, count_cut_num_of_active)
    # do loop
    kmc_loop.do_kmc_loop()
    print("###########################   normal termination   #####################################")


if __name__ == "__main__":
    main()
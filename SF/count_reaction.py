# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: count_reaction.py
@time: 12/26/2016 11:49 AM

A module to find all avail_sites for every elementary reactions.
"""


def count_of_forwards_reactions(kmc, n_avail, lat, i, j, k):
    """
    this def append '4': [['4', [3, 3], [3, 2]], ['4', [3, 3], [3, 4]]] values for keys in n_avail.

    one list element in n_avail["-0"] means one avail site for this reaction.
    The elements also contains reaction and site information.

    '4' : [['4', [j, k]], ['4', [j, k-1]]] means the reaction '4' happens on site [j, k] and [j, k-1], which means
    one site reaction.
    ['4', [3, 3], [3, 2]] means two site reaction.
    And the i store the reaction information as order in kmc.reactions

    :param kmc: instance kmc include all input information
    :param n_avail: {'4': [['4', [3, 3], [3, 2]], ['4', [3, 3], [3, 4]]]...}
    :param lat:    [['24', '24', '24', '24', '24', '24'],
                    ['18', '18', '18', '18', '18', '18'],
                    ['24', '24', '24', '24', '24', '24'],
                    ['18', '18', '18', '1',  '18', '18'],
                    ['24', '24', '24', '24', '24', '24'],
                    ['18', '18', '18', '18', '18', '18']]
    :param i: loop of all reactions
    :param j, k: loop of all surface site.
    :return: n_avail
    """
    # for case: 7 -->
    if len(kmc.reactions[i][0][0]) == 1 and kmc.reactions[i][0][0][-1].strip() in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][0][0].strip():  # if the lat[j][k] has a "7", add site info to n_avail.
            n_avail[str(i)].append([str(i), [j, k]])
    # for case 7 + CO -->
    elif len(kmc.reactions[i][0][0]) == 2 and kmc.reactions[i][0][0][-1].strip() not in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][0][0].strip():  # if the lat[j][k] has a "7", add site info to n_avail.
            n_avail[str(i)].append([str(i), [j, k]])
    # for case 7 + 8 -->
    # if the lat[j][k] has a "7", and adjacent set [j-1,k]/[j+1,k]/[j,k-1]/[j,k+1] has "8", add site info to n_avail.
    elif len(kmc.reactions[i][0][0]) == 2 and kmc.reactions[i][0][0][-1].strip() in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j-1][k] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j-1, k]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j+1][k] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k-1] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k+1] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k+1]])
    # for case 7 + 8 + CO -->
    # if the lat[j][k] has a "7", and adjacent set [j-1,k]/[j+1,k]/[j,k-1]/[j,k+1] has "8", add site info to n_avail.
    elif len(kmc.reactions[i][0][0]) == 3 and kmc.reactions[i][0][0][-1].strip() not in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j-1][k] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j-1, k]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j+1][k] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k-1] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k+1] == kmc.reactions[i][0][0][1].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k+1]])
    # for case 7 + 8 + 9 -->
    # if the lat[j][k] has a "7", and adjacent set [j-1,k]/[j+1,k] has "8"/"9", add site info to n_avail.
    # 6 cases
    elif len(kmc.reactions[i][0][0]) == 3 and kmc.reactions[i][0][0][-1].strip() in kmc.kinds:
        # for site j-1, jk, j+1,
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j-1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j-1, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j+1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j+1, k], [j-1, k]])
        # for site k-1, jk, k+1
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k-1] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j][k+1] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k-1], [j, k+1]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k+1] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j][k-1] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k+1], [j, k-1]])
        # for site j-1, jk, k-1
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j-1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j][k-1] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j-1, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k-1] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k-1], [j-1, k]])
        # for site j+1, jk, k-1
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j+1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j][k-1] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j+1, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k-1] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k-1], [j+1, k]])
        # for site j-1, jk, k+1
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j-1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j][k+1] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j-1, k], [j, k+1]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k+1] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k+1], [j-1, k]])
        # for site j+1, jk, k+1
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j+1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j][k+1] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j+1, k], [j, k+1]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j][k+1] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j, k+1], [j+1, k]])
    else:
        pass
        # print("no reaction")
    return n_avail


def count_of_reverse_reactions(kmc, n_avail, lat, i, j, k):
    """
    see notes in def count_of_forwards_reactions
    """
    # for case: <-- 7
    if len(kmc.reactions[i][0][-1]) == 1 and kmc.reactions[i][0][-1][-1].strip() in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip():  # if the lat[j][k] has a "7", add site info to n_avail.
            n_avail['-'+str(i)].append(['-'+str(i), [j, k]])
    # for case: <-- 7 + CO
    elif len(kmc.reactions[i][0][-1]) == 2 and kmc.reactions[i][0][-1][-1].strip() not in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip():  # if the lat[j][k] has a "7", add site info to n_avail.
            n_avail['-'+str(i)].append(['-'+str(i), [j, k]])
    # for case: <-- 7 + 8
    # if the lat[j][k] has a "7", and adjacent set [j-1,k]/[j+1,k]/[j,k-1]/[j,k+1] has "8", add site info to n_avail.
    elif len(kmc.reactions[i][0][-1]) == 2 and kmc.reactions[i][0][-1][-1].strip() in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k-1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k+1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k+1]])
    # for case: <-- 7 + 8 + CO
    # if the lat[j][k] has a "7", and adjacent set [j-1,k]/[j+1,k]/[j,k-1]/[j,k+1] has "8", add site info to n_avail.
    elif len(kmc.reactions[i][0][-1]) == 3 and kmc.reactions[i][0][-1][-1].strip() not in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k-1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k+1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k+1]])
    # for case: <-- 7 + 8 + 9
    # if the lat[j][k] has a "7", and adjacent set [j-1,k]/[j+1,k] has "8"/"9", add site info to n_avail.
    # 6 cases
    elif len(kmc.reactions[i][0][-1]) == 3 and kmc.reactions[i][0][-1][-1].strip() in kmc.kinds:
        # for site j-1, jk, j+1
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k], [j-1, k]])
        # for site k-1, jk, k+1
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k-1] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j][k+1] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k-1], [j, k+1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k+1] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j][k-1] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k+1], [j, k-1]])
        # for site j-1, jk, k-1
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j][k-1] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k-1] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k-1], [j-1, k]])
        # for site j+1, jk, k-1
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j][k-1] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k-1] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k-1], [j+1, k]])
        # for site j-1, jk, k+1
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j][k+1] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k], [j, k+1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k+1] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k+1], [j-1, k]])
        # for site j+1, jk, k+1
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j][k+1] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k], [j, k+1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k+1] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k+1], [j+1, k]])
    else:
        # print("no reaction")
        pass
    return n_avail


def initialize_num_of_avail_sites(kmc, lat):
    """
    def to initialize number of available sites for specific reactions.
    :param kmc: instance kmc include all input information
    :param lat:    [['24', '24', '24', '24', '24', '24'],
                    ['18', '18', '18', '18', '18', '18'],
                    ['24', '24', '24', '24', '24', '24'],
                    ['18', '18', '18', '1',  '18', '18'],
                    ['24', '24', '24', '24', '24', '24'],
                    ['18', '18', '18', '18', '18', '18']]
    :return: n_avail dict like {'22': [], '-19': [], '5': [],
                     '21': [['21', [1, 1]], ['21', [1, 2]], ['21', [1, 3]], ['21', [1, 4]], ['21', [3, 1]],
                            ['21', [3, 2]], ['21', [3, 4]]],
                     '8': [], '-23': [], '12': [], '13': [], '-5': [], '10': [], '4': [], '-14': [], '6': [],
                     '-21': [], '-24': [], '15': [], '-9': [], '-16': [], '-10': [], '11': [], '-2': [],
                     '0': [['0', [3, 3]]], '14': [], '-6': [], '-0': [], '-22': [], '7': [], '-20': [['-20', [3, 3]]],
                     '-4': [], '18': [], '-13': [], '9': [], '-15': [], '23': [], '-3': [], '20': [], '-17': [],
                     '-7': [], '3': [], '1': [['1', [3, 3]]], '2': [['2', [3, 3]]], '24': [], '-12': [], '-1': [],
                     '-18': [], '16': [], '-11': [], '17': [], '19': [], '-8': []}
    """
    n_avail = {}                         # dict to record which react happen on which site
    for i in range(len(kmc.reactions)):  # loop all reactions
        n_avail[str(i)] = []             # creat emtry list for every forward reactions
        n_avail['-'+str(i)] = []         # creat emtry list for every reverse reactions
    for i in range(len(kmc.reactions)):  # loop all reactions
        for j in range(1, len(lat)-1):   # loop all lattice, cutoff first and last one because of Periodic
            for k in range(1, len(lat)-1):  # loop all lattice, two dimension
                n_avail = count_of_forwards_reactions(kmc, n_avail, lat, i, j, k)  # add site information to reaction i
                n_avail = count_of_reverse_reactions(kmc, n_avail, lat, i, j, k)   # add site information to reaction i
    return n_avail


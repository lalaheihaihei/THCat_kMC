# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: count_reaction.py
@time: 12/26/2016 11:49 AM
"""



def count_of_forwards_reactions(kmc, n_avail, lat, i, j, k):
    """
    this def append ['-6', [3, 3], [3, 4]] values for keys in n_avail.

    one list element in n_avail["-0"] means one avail site for this reaction. the element is also an element that
    contains reaction and site information.

    "0" : [[i, j, j], [i, j, j-1]] means the reaction "0" happens on site [j, j] and [j, j-1], which means
    one site reaction and two sites reaction, respectively.
    And the i store the reaction information as order in kmc.reactions

    :param kmc: instance kmc include all input information
    :param n_avail: {'-6': [['-6', 3, 3]], '5': [], '-2': []...}
    :param lat: [['0', '0', '0', '1', '0', '0'],['0','0','0','0','0','0']]
    :param i: loop of all reactions
    :param j: loop of all surface site. In oneD module, only one site or two site reaction are allowed.
    :return: n_avail
    """
    # for case: 7 -->
    if len(kmc.reactions[i][0][0]) == 1 and kmc.reactions[i][0][0][0].strip() in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][0][0].strip():
            n_avail[str(i)].append([str(i), [j, k]])
    # for case 7 + CO -->
    elif len(kmc.reactions[i][0][0]) == 2 and kmc.reactions[i][0][0][-1].strip() not in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][0][0].strip():
            n_avail[str(i)].append([str(i), [j, k]])
    # for case 7 + 8 -->
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
    elif len(kmc.reactions[i][0][0]) == 3 and kmc.reactions[i][0][0][-1].strip() in kmc.kinds:
        # for site j-1, j, j+1
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j-1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j-1, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][0][0].strip() and lat[j+1][k] == kmc.reactions[i][0][0][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][0][2].strip():
            n_avail[str(i)].append([str(i), [j, k], [j+1, k], [j-1, k]])
        # for site k-1, k, k+1
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
    # for case: 7 -->
    if len(kmc.reactions[i][0][-1]) == 1 and kmc.reactions[i][0][-1][0].strip() in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k]])
    # for case 7 + CO -->
    elif len(kmc.reactions[i][0][-1]) == 2 and kmc.reactions[i][0][-1][-1].strip() not in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k]])
    # for case 7 + 8 -->
    elif len(kmc.reactions[i][0][-1]) == 2 and kmc.reactions[i][0][-1][-1].strip() in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k-1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k+1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k+1]])
    # for case 7 + 8 + CO -->
    elif len(kmc.reactions[i][0][-1]) == 3 and kmc.reactions[i][0][-1][-1].strip() not in kmc.kinds:
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k-1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k-1]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j][k+1] == kmc.reactions[i][0][-1][1].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j, k+1]])
    # for case 7 + 8 + 9 -->
    elif len(kmc.reactions[i][0][-1]) == 3 and kmc.reactions[i][0][-1][-1].strip() in kmc.kinds:
        # for site j-1, j, j+1
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j-1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j+1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j-1, k], [j+1, k]])
        if lat[j][k] == kmc.reactions[i][0][-1][0].strip() and lat[j+1][k] == kmc.reactions[i][0][-1][1].strip() and \
                        lat[j-1][k] == kmc.reactions[i][0][-1][2].strip():
            n_avail['-'+str(i)].append(['-'+str(i), [j, k], [j+1, k], [j-1, k]])
        # for site k-1, k, k+1
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
    :param lat: [['0', '0', '0', '1', '0', '0'],['0','0','0','0','0','0']]
    :return: n_avail dict like {'-6': [['-6', 3, 3]], '5': [], '-2': []...}
    """
    n_avail = {}                        # dict to record which react happen on which site {'-0':[1,5,8,...],'0':[]}
    for i in range(len(kmc.reactions)):
        n_avail[str(i)] = []            # creat list to record site
        n_avail['-'+str(i)] = []        # creat list to record site
    for i in range(len(kmc.reactions)): # loop all reactions
        for j in range(1, len(lat)-1):  # loop all lattice ,cutoff first and last to avoid double count in Periodic case
            for k in range(1, len(lat)-1):
                n_avail = count_of_forwards_reactions(kmc, n_avail, lat, i, j, k)
                n_avail = count_of_reverse_reactions(kmc, n_avail, lat, i, j, k)
    return n_avail


# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py
@time: 11/29/2016 11:15 AM
main program to run kMC
"""
import math
import random
import sys

from SAmodule import adsorption
from SAmodule import reactions
from SAmodule import species
from IO import Config
from SAmodule import lattice


def test_kinds_number(kmc):
    """
    Test the input species number. Whether the TS + Local minimums = tot species.

    :param kmc: Instance of Parameters
    :return: None
    """
    if kmc.isKindsTS.count(0) != len(kmc.kinds) - kmc.isKindsTS.count(1):
        raise ValueError('Error, please checks the kinds and isKindsTS parameters.')
    return None


def test_reaction_number(kmc):
    """
    Test the input reaction number. Whether the n(ads+des+rea) = n(reactions)

    :param kmc: Instance of Parameters
    :return: None
    """
    if len(kmc.reactionsKind) != len(kmc.reactions):
        print(kmc.reactionsKind, kmc.reactions)
        raise ValueError('Error: Please check the reactions number and reactionsKind number')
    return None


def add_list(kmc, speciesList = [], TSList = [], totList = []):
    """
    Set species class with kind, frequency, and energy information
    and append them to speciesList, TSList, and totList.

    :param kmc: Instance of Parameters
    :param speciesList: New List to store the local minimum species information.
    :param TSList: New List to store the TS species information.
    :param totList: New List to store the total species information.
    :return speciesTuple: Tuple to store the local minimum species information.
    :return speciesTuple_kind: Tuple to store the local minimum species kinds.
    :return TSTuple: Tuple to store the TS species information.
    :return totList: Tuple to store the total species information.
    :return totList_kind: Tuple to store the total species kinds.

    """
    for i in range(len(kmc.kinds)):
        if kmc.isKindsTS[i] == 0:
            #print(kmc.kinds, kmc.kindsFreq, kmc.kindsEnergy)
            speciesList.append(species.species(kmc.kinds[i], kmc.kindsFreq[i], kmc.kindsEnergy[i]))
            totList.append(speciesList[-1])
        elif kmc.isKindsTS[i] == 1:
            TSList.append(species.species(kmc.kinds[i], kmc.kindsFreq[i], kmc.kindsEnergy[i]))
            totList.append(TSList[-1])
        else:
            raise ValueError('Error: Error, please checks the isKindsTS parameters.')
    speciesTuple = tuple(speciesList)
    speciesTuple_kind = tuple(i.kind for i in speciesTuple)
    TSTuple = tuple(TSList)
    totList_kind = tuple(i.kind for i in totList)
    return speciesTuple, speciesTuple_kind, TSTuple, totList, totList_kind


def add_reaction(kmc, totList, T, P, reactionList = [], k_forward = [], k_reverse = [],\
                 dict_reaction = {}, dict_reaction_direction = {}):
    """
    Calculation the reaction rate and store reaction informations.

    :param kmc: Instance of Parameters
    :param totList: Tuple with the total species information.
    :param T: Temperature
    :param P: Partial pressure
    :param reactionList: New List to store all the reaction informations.
    :param k_forward: New List to store all the forward reaction k.
    :param k_reverse: New List to store all the reverse reaction k.
    :param dict_reaction: New dict to store {IS:[k]}, ex:
        {0: [12.934191, 12.93419, 10.318585], 1: [0.0004196, 0.04606873]}
    :param dict_reaction_direction: New dict to store {IS:[FS]}, ex:
        {0: [1, 9, 8], 1: [0, 3], 3: [1, 4], 4: [3, 5, 11], 5: [4, 6]}
    :return: reactionList, k_forward, k_reverse, dict_reaction, dict_reaction_direction
    """
    for i in range(len(kmc.reactions)):
        # judge the reaction kind, and calculate k and G, append it to reactionList.
        if kmc.reactionsKind[i] == 'ads':
            reactionList.append(adsorption.Adsorption\
                (totList[int(kmc.reactions[i][0][0])], totList[int(kmc.reactions[i][1][0])], T, \
                 P[kmc.reactions[i][0][1].strip()], kmc.reactions[i][0][1].strip()))
            print("For step %d:\tk(ads) = %.3e,\tk(des) = %.3e,\tDeltaG = %.3f" % \
                  (i + 1, reactionList[-1].adsorb_k(), reactionList[-1].desorb_k(), reactionList[-1].delta_g()))
            k_forward.append(reactionList[-1].adsorb_k())  # WHY does work by directly append like k_forward.append(reactionList[-1].adsorbK())
            k_reverse.append(reactionList[-1].desorb_k())

            # if IS in dict_reaction_direction: append(k), if not create IS:[k]
            if int(kmc.reactions[i][0][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][0][0])].append(k_forward[-1])
            else:
                dict_reaction[int(kmc.reactions[i][0][0])] = [k_forward[-1]]
            if int(kmc.reactions[i][1][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][1][0])].append(k_reverse[-1])
            else:
                dict_reaction[int(kmc.reactions[i][1][0])] = [k_reverse[-1]]

            # if IS in dict_reaction_direction: append(FS), if not create IS:[FS]
            if int(kmc.reactions[i][0][0]) in dict_reaction_direction.keys():
                dict_reaction_direction[int(kmc.reactions[i][0][0])].append(int(kmc.reactions[i][1][0]))
            else:
                dict_reaction_direction[int(kmc.reactions[i][0][0])] = [int(kmc.reactions[i][1][0])]
            if int(kmc.reactions[i][1][0]) in dict_reaction_direction.keys():
                dict_reaction_direction[int(kmc.reactions[i][1][0])].append(int(kmc.reactions[i][0][0]))
            else:
                dict_reaction_direction[int(kmc.reactions[i][1][0])] = [int(kmc.reactions[i][0][0])]

        elif kmc.reactionsKind[i] == 'react':
            reactionList.append(reactions.reactions\
                    (totList[int(kmc.reactions[i][0][0])], totList[int(kmc.reactions[i][1][0])], \
                     totList[int(kmc.reactions[i][2][0])], T))
            print("For step %d:\tk(forwards) = %.3e,\tk(reverse) = %.3e"\
                  % (i+1, reactionList[-1].forwardK, reactionList[-1].reverseK))
            k_forward.append(reactionList[-1].forwardK)
            k_reverse.append(reactionList[-1].reverseK)
            if int(kmc.reactions[i][0][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][0][0])].append(k_forward[-1])
            else:
                dict_reaction[int(kmc.reactions[i][0][0])] = [k_forward[-1]]
            if int(kmc.reactions[i][2][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][2][0])].append(k_reverse[-1])
            else:
                dict_reaction[int(kmc.reactions[i][2][0])] = [k_reverse[-1]]

            if int(kmc.reactions[i][0][0]) in dict_reaction_direction.keys():
                dict_reaction_direction[int(kmc.reactions[i][0][0])].append(int(kmc.reactions[i][2][0]))
            else:
                dict_reaction_direction[int(kmc.reactions[i][0][0])] = [int(kmc.reactions[i][2][0])]
            if int(kmc.reactions[i][2][0]) in dict_reaction_direction.keys():
                dict_reaction_direction[int(kmc.reactions[i][2][0])].append(int(kmc.reactions[i][0][0]))
            else:
                dict_reaction_direction[int(kmc.reactions[i][2][0])] = [int(kmc.reactions[i][0][0])]

        elif kmc.reactionsKind[i] == 'des':
            reactionList.append(adsorption.Adsorption\
                (totList[int(kmc.reactions[i][1][0])], totList[int(kmc.reactions[i][0][0])], T, \
                 P[kmc.reactions[i][1][1].strip()], kmc.reactions[i][1][1].strip()))
            print("For step %d:\tk(des) = %.3e,\tk(ads) = %.3e,\tDeltaG = %.3f" % \
                  (i + 1, reactionList[-1].desorb_k(), reactionList[-1].adsorb_k(), - reactionList[-1].delta_g()))
            k_forward.append(reactionList[-1].desorb_k())
            k_reverse.append(reactionList[-1].adsorb_k())
            if int(kmc.reactions[i][0][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][0][0])].append(k_forward[-1])
            else:
                dict_reaction[int(kmc.reactions[i][0][0])] = [k_forward[-1]]
            if int(kmc.reactions[i][1][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][1][0])].append(k_reverse[-1])
            else:
                dict_reaction[int(kmc.reactions[i][1][0])] = [k_reverse[-1]]
            #print(reactionList[-1].gas_mu())
            if int(kmc.reactions[i][0][0]) in dict_reaction_direction.keys():
                dict_reaction_direction[int(kmc.reactions[i][0][0])].append(int(kmc.reactions[i][1][0]))
            else:
                dict_reaction_direction[int(kmc.reactions[i][0][0])] = [int(kmc.reactions[i][1][0])]
            if int(kmc.reactions[i][1][0]) in dict_reaction_direction.keys():
                dict_reaction_direction[int(kmc.reactions[i][1][0])].append(int(kmc.reactions[i][0][0]))
            else:
                dict_reaction_direction[int(kmc.reactions[i][1][0])] = [int(kmc.reactions[i][0][0])]


        else:
            raise ValueError('Error: please check the reactions parameters.')
    # Add the last step to the first one: Warming, the last reaction counted twice.

    dict_reaction[0].append(dict_reaction[len(totList)-1][0])
    dict_reaction_direction[0].append(dict_reaction_direction[len(totList)-1][0])
    return reactionList, k_forward, k_reverse, dict_reaction, dict_reaction_direction


def initialize_lattice(kmc, speciesTuple):
    """
    initialize the lattice information: Maybe no use for SA simualtions.
    :param kmc: Instance of Parameters
    :param speciesTuple: Tuple to store the local minimum species information.
    :return:
    """
    lat = lattice.lattice(kmc.latticeSize)    # set lattice dimension.
    print('lattice size is setted to be: %d * %d' %(lat.demension[0], lat.demension[1]))
    lat.initialElements = speciesTuple[0].kind    # choose the initial kind on surface.
    # print('lattice elements is initialized to be %s' % (lat.initialElements))
    lat.initializeElements()    # initialize the elements on surface.
    print('lattice is initialized')
    return lat


def do_kmc_loop(kmc, lat, speciesTuple_kind, totList, totList_kind, dict_reaction, dict_reaction_direction,\
        time = [0], productNum = 0, coverage_t_List = []):
    """
    run kMC simulation.

    :param kmc: Instance of Parameters
    :param lat: lattice information
    :param speciesTuple_kind: Tuple to store the local minimum species kinds.
    :param totList: Tuple to store the total species information.
    :param totList_kind: Tuple to store the total species kinds.
    :param dict_reaction: Dict to store {IS:[k]}
    :param dict_reaction_direction:  Dict to store {IS:[FS]}
    :param time: initilize the start time to 0
    :param productNum: count of the products
    :param coverage_t_List: Add a time dependent coverage
    :return: time, productNum, coverageList
    """
    for i in range(len(speciesTuple_kind)):
        coverage_t_List.append(0)
    for i in range(kmc.LoopNum):                       # MC cycle number
        k_rate = k_rate1 = 0                           # initialize k
        theReactionSite = []
        for j in range(lat.demension[0]):              # x direction
            for k in range(lat.demension[1]):          # y direction
                IS = lat.elements[k][j]                # find out all the IS on surface
                IS_order = totList_kind.index(IS)    # find the IS's corresponding IS_order
                if IS_order == len(kmc.kinds):                         # s(last) = s1
                    IS_order = 0
                #print(IS_order, IS_order-1, k_forward[IS_order], k_reverse[IS_order-1])
                k_rate += sum(dict_reaction[IS_order])  # get K_tot


        # bug fixed, coverage information should be recorded firstly. Then, update the reaction directon.
        # update time
        t_random = random.random()
        d_t = (1/k_rate) * math.log(1/t_random)
        time.append(time[-1] + d_t)

        #add a time dependent coverage: see Top Catal (2014) 57:159–170 equation (6).
        coverage_t_List = list(map( lambda i : coverage_t_List[i] + d_t * lat.elements[0].count(speciesTuple_kind[i])/\
                lat.demension[0] * lat.demension[1], [ i for i in range(len(speciesTuple_kind))]))



        # bug fixed, coverage information should be recorded firstly. Then, update the reaction directon.
        k_random = random.random()                 # random a k_random
        # print('k_random and k_rate =', k_random,k_rate)
        for j in range(lat.demension[0]):       # x direction
            for k in range(lat.demension[1]):   # y direction
                IS = lat.elements[k][j]         # find out all the IS on surface
                IS_order = totList_kind.index(IS)
                if IS_order == len(totList) - 1:
                    IS_order = 0
                    # if IS is the last one in speciesTuple'\
                    # make the it to be the first one in speciesTupl.
                for l in range(len(dict_reaction[IS_order])):
                    k_rate1 += dict_reaction[IS_order][l]
                    if k_rate1 > k_random * k_rate:                                  # get the reaction site
                        theReactionSite.append(j)
                        theReactionSite.append(k)
                        theReactionSite.append(dict_reaction_direction[IS_order][l])          # get the reaction direction
                        #print("break")
                        break
                else:
                    continue
                break
            else:
                continue
            break
        #print(theReactionSite)

        FS_order = theReactionSite[2]   # update FS_order


        '''
        if FS_order == -1:
            FS_order = len(speciesTuple) - 2
            # correct the FS_order, if FS_order = -1. \
            # which means IS = speciesTuple[0] --> FS = speciesTuple[-2]
        '''
        if FS_order == len(totList) - 1:
            FS_order = 0
            # correct the FS_order, if FS_order = len(speciesTuple) - 1. \
            # which means IS = speciesTuple[6] --> FS = speciesTuple[7]
        #print(IS_order, FS_order)

        FS = totList_kind[FS_order]

        # count the reaction sites:
        if totList[FS_order] == totList[0] and totList[IS_order] == totList[-2]:
            #print('********* +1')
            productNum += 1
        elif totList[FS_order] == totList[-2] and totList[IS_order] == totList[0]:
            #print('********* -1')
            productNum -= 1
        else:
            pass

        lat.elements = ((theReactionSite[0], theReactionSite[1]), IS, FS) # update lattice elements


    coverageList = list(map(lambda i: i / time[-1], coverage_t_List))

    return time, productNum, coverageList


def main():

    kmc = Config.Parameters()

    T = kmc.Temperature                                               # get temperature /K
    P = {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))}  # get pressure    /Bar

    # Test the input kinds and reaction number.
    test_kinds_number(kmc)
    test_reaction_number(kmc)

    # Input all species information.
    speciesTuple, speciesTuple_kind, TSTuple, totList, totList_kind = add_list(kmc)

    # Input all reactions information.
    reactionList, k_forward, k_reverse, dict_reaction, dict_reaction_direction\
        = add_reaction(kmc, totList, T, P)

    # Output the reactions information.
    print("k_forward and k_reverse :", k_forward, '\n', k_reverse)
    print("dict_reaction :", dict_reaction)
    print("dict_reaction_direction :", dict_reaction_direction)

    # Initialize the lattice information. (Not necessary for SA part)
    lat = initialize_lattice(kmc, speciesTuple)

    # DO KMC
    time, productNum, coverageList = do_kmc_loop(kmc, lat, speciesTuple_kind, totList, totList_kind, dict_reaction, dict_reaction_direction)

    # print(lat.elements)

    # print coverage information, time, and product number.
    print(speciesTuple, coverageList)
    for i in range(len(speciesTuple_kind)):
        print('the coverage for species {:s} = {:.3f}'. format(speciesTuple[i].kind, coverageList[i]))
    print('kMC simulation time is %.2f' % time[-1], 'product number is %d ' % productNum)

    # Output TOF.
    print('TOF calculated by product/(time * sites) by ln(r/s-1) = %.2f ' % math.log(productNum / (time[-1] * lat.demension[0]) * int(lat.demension[1])))

    # this method may have some problem, (1-r) should be multiple as discussed in Norskov'book: (7.4)
    print('TOF calculated by Theta * k_forward by ln(r/s-1) = %.2f ' \
          % math.log(max(coverageList) * k_forward[coverageList.index(max(coverageList))] ) )

if __name__ == "__main__":
    main()
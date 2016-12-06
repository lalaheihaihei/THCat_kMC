# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py
@time: 11/29/2016 11:15 AM
main program to run kMC
"""
import math
import random

from SAmodule import adsorption
from SAmodule import reactions
from SAmodule import species
from IO import Config
from SAmodule import lattice


def test_kinds_number(kmc):
    if kmc.isKindsTS.count(0) != len(kmc.kinds) - kmc.isKindsTS.count(1):
        raise ValueError('Error, please checks the kinds and isKindsTS parameters.')
    return None


# Add set species class with kind, frequency, and energy information and append them to speciesList.
def add_list(kmc, speciesList = [], TSList = [], totList = []):
    for i in range(len(kmc.kinds)):
        if kmc.isKindsTS[i] == 0:
            speciesList.append(species.species(kmc.kinds[i], kmc.kindsFreq[i], kmc.kindsEnergy[i]))
            totList.append(speciesList[-1])
        elif kmc.isKindsTS[i] == 1:
            TSList.append(species.species(kmc.kinds[i], kmc.kindsFreq[i], kmc.kindsEnergy[i]))
            totList.append(TSList[-1])
        else:
            raise ValueError('Error: Error, please checks the isKindsTS parameters.')
    # Transfer local minimum list and ts list to tuple format.
    speciesTuple = tuple(speciesList)
    speciesTuple_kind = tuple(i.kind for i in speciesTuple)
    TSTuple = tuple(TSList)
    totList_kind = tuple(i.kind for i in totList)
    # print('totlist = ', totList)
    return speciesTuple, speciesTuple_kind, TSTuple, totList, totList_kind

# print(speciesTuple, speciesTuple_kind, TSTuple)
# speciesTuple[8] = speciesTuple[1]

def add_reaction(kmc, totList, T, P, reactionList = [], k_forward = [], k_reverse = [], dict_reaction = {}, dict_reaction_direction = {}):
    for i in range(len(kmc.reactions)):
        if kmc.reactionsKind[i] == 'ads':
            reactionList.append(adsorption.Adsorption\
                (totList[int(kmc.reactions[i][0][0])], totList[int(kmc.reactions[i][1][0])], T, \
                 P[kmc.reactions[i][0][1].strip()], kmc.reactions[i][0][1].strip()))
            print("For step %d:\tk(ads) = %.3e,\tk(des) = %.3e,\tDeltaG = %.3f" % \
                  (i + 1, reactionList[-1].adsorb_k(), reactionList[-1].desorb_k(), reactionList[-1].delta_g()))
            f = reactionList[-1].adsorb_k()
            r = reactionList[-1].desorb_k()
            k_forward.append(f)  # WHY does work by directly append like k_forward.append(reactionList[-1].adsorbK())
            k_reverse.append(r)
            if int(kmc.reactions[i][0][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][0][0])].append(k_forward[-1])
            else:
                dict_reaction[int(kmc.reactions[i][0][0])] = [k_forward[-1]]
            if int(kmc.reactions[i][1][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][1][0])].append(k_reverse[-1])
            else:
                dict_reaction[int(kmc.reactions[i][1][0])] = [k_reverse[-1]]

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
            f = reactionList[-1].forwardK
            r = reactionList[-1].reverseK
            k_forward.append(f)
            k_reverse.append(r)
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
            f = reactionList[-1].desorb_k()
            r = reactionList[-1].adsorb_k()
            k_forward.append(f)
            k_reverse.append(r)
            if int(kmc.reactions[i][0][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][0][0])].append(k_forward[-1])
            else:
                dict_reaction[int(kmc.reactions[i][0][0])] = [k_forward[-1]]
            if int(kmc.reactions[i][1][0]) in dict_reaction.keys():
                dict_reaction[int(kmc.reactions[i][1][0])].append(k_reverse[-1])
            else:
                dict_reaction[int(kmc.reactions[i][1][0])] = [k_reverse[-1]]

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
    return reactionList, k_forward, k_reverse, dict_reaction, dict_reaction_direction


def initialize_lattice(kmc, speciesTuple):
    lat = lattice.lattice(kmc.latticeSize)    # set lattice dimension.
    print('lattice size is setted to be: %d * %d' %(lat.demension[0], lat.demension[1]))
    lat.initialElements = speciesTuple[0].kind    # choose the initial kind on surface.
    print('lattice elements is initialized to be %s' % (lat.initialElements))
    lat.initializeElements()    # initialize the elements on surface.
    print('lattice is initialized')
    return lat


def main():

    kmc = Config.Parameters()

    # Basic parameters.
    T = kmc.Temperature                                               # get temperature /K
    P = {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))}  # get pressure    /Bar
    latticeSize = kmc.latticeSize


    # Input all species information.
    localMinimum = kmc.isKindsTS.count(0)
    test_kinds_number(kmc)
    speciesTuple, speciesTuple_kind, TSTuple, totList, totList_kind = add_list(kmc)



    # print(kmc.reactions)
    # find out all the pathways and rate constant.
    reactionList, k_forward, k_reverse, dict_reaction, dict_reaction_direction = add_reaction(kmc, totList, T, P)

    # record the reaction direction correspond the dict_reaction
    if len(kmc.reactionsKind) != len(kmc.reactions):
        raise ValueError('Error: Please check the reactions number and reactionsKind number')


    # reactionTuple = tuple(reactionList)
    print(reactionList)
    print(k_forward, '\n', k_reverse)
    dict_reaction[0].append(dict_reaction[len(totList)-1][0])
    dict_reaction_direction[0].append(dict_reaction_direction[len(totList)-1][0])
    print(dict_reaction)
    print(dict_reaction_direction)

    # DO lattice KMC
    lat = initialize_lattice(kmc, speciesTuple)

    time = [0]                                         # initialize time
    productNum = 0                                     # initialize the product number
    coverage_t_List = []
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
        if totList[FS_order] == totList[0] and totList[IS_order] == totList[8]:
            #print('********* +1')
            productNum += 1
        elif totList[FS_order] == totList[8] and totList[IS_order] == totList[0]:
            #print('********* -1')
            productNum -= 1
        else:
            pass

        lat.elements = ((theReactionSite[0], theReactionSite[1]), IS, FS) # update lattice elements

        # update time
        t_random = random.random()
        d_t = (1/k_rate) * math.log(1/t_random)
        time.append(time[-1] + d_t)

        #add a time dependent coverage: see Top Catal (2014) 57:159–170 equation (6).
        coverage_t_List = list(map( lambda i : coverage_t_List[i] + d_t * lat.elements[0].count(speciesTuple_kind[i])/\
                lat.demension[0] * lat.demension[1], [ i for i in range(len(speciesTuple_kind))]))

    coverageList = list(map( lambda i : i / time[-1] , coverage_t_List))

    #print(lat.elements)
    for i in range(len(speciesTuple_kind)):
        print('the coverage for species {:s} = {:.3f}'. format(speciesTuple[i].kind, coverageList[i]))
    print('kMC simulation time is %.2f' % time[-1], 'product number is %d ' % productNum)
    print('TOF calculated by product/(time * sites) by ln(r/s-1) = %.2f ' % math.log(productNum / (time[-1] * lat.demension[0]) * int(lat.demension[1])))

    # this method may have some problem, (1-r) should be multiple as discussed in Norskov'book: (7.4)
#    print('TOF calculated by Theta * k_forward by ln(r/s-1) = %.2f ' \
#          % math.log(max(coverageList) * k_forward[coverageList.index(max(coverageList))] ) )

if __name__ == "__main__":
    main()
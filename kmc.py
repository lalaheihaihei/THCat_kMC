# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py
@time: 11/29/2016 11:15 AM
main program to run kMC
"""
import species, lattice, reactions, adsorption, re, math, random, yaml, configparser, Config
import scipy.constants as sc

kmc = Config.Parameters()

# Basic parameters.
T = kmc.Temperature                                               # get temperature /K
P = {kmc.gas[i]:kmc.GasPressure[i] for i in range(len(kmc.gas))}  # get pressure    /Bar
latticeSize = kmc.latticeSize



# Input all species information.
localMinimum = kmc.isKindsTS.count(0)
speciesList = []
TSList = []
totList = []
if localMinimum != len(kmc.kinds) - kmc.isKindsTS.count(1):
    raise ValueError('Error, please checks the kinds and isKindsTS parameters.')
# Add set species class with kind, frequency, and energy information and append them to speciesList.
for i in range(len(kmc.kinds)):
    if kmc.isKindsTS[i] == 0:
        speciesList.append(species.species(kmc.kinds[i], kmc.kindsFreq[i], kmc.kindsEnergy[i]))
        totList.append(speciesList[-1])
    elif kmc.isKindsTS[i] == 1:
        TSList.append(species.species(kmc.kinds[i], kmc.kindsFreq[i],kmc.kindsEnergy[i]))
        totList.append(TSList[-1])
    else:
        raise ValueError('Error: Error, please checks the isKindsTS parameters.')
# Transfer local minimum list and ts list to tuple format.
speciesTuple = tuple(speciesList)
speciesTuple_kind = tuple(i.kind for i in speciesTuple)
TSTuple = tuple(TSList)
totList = tuple(totList)
#print(speciesTuple, speciesTuple_kind, TSTuple)         #speciesTuple[8] = speciesTuple[1]
# Abbreviate speciesTuple to s, and TSTuple to ts.
s = totList


#print(kmc.reactions)
# find out all the pathways and rate constant.

reactionList = []
k_forward = []
k_reverse = []
if len(kmc.reactionsKind) != len(kmc.reactions):
    raise ValueError('Error: Please check the reactions number and reactionsKind number')

for i in range(len(kmc.reactions)):
    if kmc.reactionsKind[i] == 'ads':
        reactionList.append(adsorption.adsorption\
            (s[int(kmc.reactions[i][0][0])], s[int(kmc.reactions[i][1][0])], T,\
             P[kmc.reactions[i][0][1].strip()], kmc.reactions[i][0][1].strip()))
        print("For step %d:\tk(ads) = %.3e,\tk(des) = %.3e,\tDeltaG = %.3f" % \
              (i+1, reactionList[-1].adsorbK(), reactionList[-1].desorbK(), reactionList[-1].deltaG()))
        f = reactionList[-1].adsorbK()
        r = reactionList[-1].desorbK()
        k_forward.append(f)  # WHY does work by directly append like k_forward.append(reactionList[-1].adsorbK())
        k_reverse.append(r)
    elif kmc.reactionsKind[i] == 'react':
        reactionList.append(reactions.reactions\
                (s[int(kmc.reactions[i][0][0])], s[int(kmc.reactions[i][1][0])],\
                 s[int(kmc.reactions[i][2][0])], T))
        print("For step %d:\tk(forwards) = %.3e,\tk(reverse) = %.3e"\
              % (i+1, reactionList[-1].forwardK, reactionList[-1].reverseK))
        f = reactionList[-1].forwardK
        r = reactionList[-1].reverseK
        k_forward.append(f)
        k_reverse.append(r)
    elif kmc.reactionsKind[i] == 'des':
        reactionList.append(adsorption.adsorption\
            (s[int(kmc.reactions[i][1][0])], s[int(kmc.reactions[i][0][0])], T,\
             P[kmc.reactions[i][1][1].strip()], kmc.reactions[i][1][1].strip()))
        print("For step %d:\tk(des) = %.3e,\tk(ads) = %.3e,\tDeltaG = %.3f" % \
              (i+1, reactionList[-1].desorbK(), reactionList[-1].adsorbK(), reactionList[-1].deltaG()))
        f = reactionList[-1].desorbK()
        r = reactionList[-1].adsorbK()
        k_forward.append(f)
        k_reverse.append(r)
    else:
        raise ValueError('Error: please check the reactions parameters.')
#reactionTuple = tuple(reactionList)
print(reactionList)
print(k_forward, '\n', k_reverse)

# DO lattice KMC
lat = lattice.lattice(kmc.latticeSize)    # set lattice dimension.
print('lattice size is setted to be: %d * %d' %(lat.demension[0], lat.demension[1]))
lat.initialElements = speciesTuple[0].kind    # choose the initial kind on surface.
print('lattice elements is initialized to be %s' % (lat.initialElements))
lat.initializeElements()    # initialize the elements on surface.
print('lattice is initialized')

time = [0]                                           # initialize time
productNum = 0                                      # initialize the product number
for i in range(kmc.LoopNum):                       # MC cycle number
    k_rate = k_rate1 = 0                           # initialize k
    theReactionSite = []
    for j in range(lat.demension[0]):              # x direction
        for k in range(lat.demension[1]):          # y direction
            IS = lat.elements[k][j]                # find out all the IS on surface
            IS_order = speciesTuple_kind.index(IS)    # find the IS's corresponding IS_order
            if IS_order == len(kmc.kinds):                         # s8 = s1
                IS_order = 0
            #print(IS_order, IS_order-1, k_forward[IS_order], k_reverse[IS_order-1])
            k_rate += k_forward[IS_order] + k_reverse[IS_order-1]  # get K_tot

    k_random = random.random()                 # random a k_random
    #print(k_random,k_rate)
    for j in range(lat.demension[0]):       # x direction
        for k in range(lat.demension[1]):   # y direction
            IS = lat.elements[k][j]         # find out all the IS on surface
            IS_order = speciesTuple_kind.index(IS)
            if IS_order == len(speciesTuple) - 1:
                IS_order = 0
                # if IS is the last one in speciesTuple'\
                # make the it to be the first one in speciesTupl.
            k_rate1 += k_forward[IS_order]
            if  k_rate1 < k_random * k_rate:                  # once random > tot,
                pass
            elif k_rate1 > k_random * k_rate:                                  # get the reaction site
                theReactionSite.append(j)
                theReactionSite.append(k)
                theReactionSite.append(1)          # get the reaction direction
                #print("break")
                break

            k_rate1 += k_reverse[IS_order - 1]
            if  k_rate1 < k_random * k_rate:
                pass
            elif k_rate1 > k_random * k_rate:
                theReactionSite.append(j)
                theReactionSite.append(k)
                theReactionSite.append(-1)
                #print("break")
                break
        else:
            continue
        break
    #print(theReactionSite)
    FS_order = IS_order + theReactionSite[2]   # update FS_order

    if FS_order == -1:
        FS_order = len(speciesTuple) - 2
        # correct the FS_order, if FS_order = -1. \
        # which means IS = speciesTuple[0] --> FS = speciesTuple[-2]
    elif FS_order == len(speciesTuple) - 1:
        FS_order = 0
        # correct the FS_order, if FS_order = len(speciesTuple) - 1. \
        # which means IS = speciesTuple[6] --> FS = speciesTuple[7]
    #print(IS_order, FS_order)

    FS = speciesTuple[FS_order].kind

    # count the reaction sites:
    if speciesTuple[FS_order] == speciesTuple[0] and speciesTuple[IS_order] == speciesTuple[-2]:
        #print('**** +1')
        productNum += 1
    elif speciesTuple[FS_order] == speciesTuple[-2] and speciesTuple[IS_order] == speciesTuple[0]:
        #print('**** -1')
        productNum -= 1
    else:
        pass

    lat.elements = ((theReactionSite[0], theReactionSite[1]), IS, FS) # update lattice elements

    # update time
    t_random = random.random()
    d_t = (1/k_rate) * math.log(1/t_random)
    time.append(time[-1] + d_t)

#print(lat.elements)
for i in range(len(speciesTuple_kind)):
    print('the coverage for species {:d} = {:.3f}'.\
          format(i, float(lat.elements[0].count(speciesTuple_kind[i])) / \
            int(lat.demension[0]) * int(lat.demension[1])))
print('kMC simulation time is %.2f' % time[-1], 'product number is %d ' % productNum)


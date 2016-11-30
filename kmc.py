# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc.py
@time: 11/29/2016 11:15 AM
main program to run kMC
"""
import species, lattice, reactions, adsorption, re, math, random, yaml, configparser, Config
import scipy.constants as sc

'''
TEST CODES:

co = species.species('co',(123,2000),1000)

print("%2.3f" % (co.ZPE))

CeO2_stepD = lattice.lattice((10,10))

print(CeO2_stepD.demension)

CeO2_stepD.initialElements = '*'

print(CeO2_stepD.initialElements)

CeO2_stepD.initializeElements()

site = (8,0)
IS = '*'
FS = 'CO'
print('the \'%s\' on %d,%d site was replaced by \'%s\'' % (IS, site[0], site[1], FS))
CeO2_stepD.elements = ((8,0),'*','CO')
print(CeO2_stepD.elements)



co = species.species('co',(123,2000),1000)

print("%2.3f" % (co.correctEnergy))

print(sc.physical_constants['Boltzmann constant in eV/K'])
print(sc.physical_constants['Planck constant in eV s'])


IS = species.species('co',(123,2000),999)
TS = species.species('co-ts-co2',(123,2000),1000)
FS = species.species('co2',(123,2000),998)
T  = 300
path1 = reactions.reactions(IS, TS, FS, T)
print(path1.forwardK,path1.reverseK)
'''

'''
OLD TEST CODES for kMC, just keep here for comparation:

# parameters
T = 300 # temperature K
P_CO = 1e-4 # paritial pressure of CO
P_O2 = 1e-2
P_CO2 = 1e-4
latticeSize = (100,1)

# find all the reaction pathway
s1 = species.species('*OAu    + *O',(123,2000),0.00)
s2 = species.species('*OAu-CO + *O',(123,2000),-1.308)
ts2_3 = species.species('*OAu-CO...O*',(123,2000),-0.479)
s3 = species.species('*OAu-COO*',(123,2000),-1.240)
s4 = species.species('*OAu    + *',(123,2000),-1.140)
s5 = species.species('*OAu    + *O2',(123,2000),-2.709)
s6 = species.species('*OAuCO  + *O2',(123,2000),-4.155)
ts6_7 = species.species('*OAuCO...O2*',(123,2000),-3.391)
s7 = species.species('*OAu-CO-OO*',(123,2000),-5.412)
s8 = species.species('*OAu    + *O',(123,2000),-6.544)
#s8 = s1
speciesTuple = (s1,s2,s3,s4,s5,s6,s7,s8)
speciesTuple_kind = tuple(i.kind for i in speciesTuple)

# find out rate constant
path1_2 = adsorption.adsorption(s1, s2, T, P_CO, 'CO')
print("For step one:\tk(ads) = %.3e,\tk(des) = %.3e,\tDeltaG = %.3f\n" \
      % (path1_2.adsorbK(), path1_2.desorbK(), path1_2.deltaG()))

path2_3 = reactions.reactions(s2, ts2_3, s3, T)
print("For step two:\tk(forwards) = %.3e,\tk(reverse) = %.3e,\n"%(path2_3.forwardK,path2_3.reverseK))

path3_4 = adsorption.adsorption(s4, s3, T, P_CO2, 'CO2')
print("For step three:\tk(des) = %.3e,\tk(ads) = %.3e,\tDeltaG = %.3f\n" \
      % (path3_4.desorbK(), path3_4.adsorbK(), - path3_4.deltaG()))

path4_5 = adsorption.adsorption(s4, s5, T, P_O2, 'O2')
print("For step four:\tk(ads) = %.3e,\tk(des) = %.3e,\tDeltaG = %.3f\n" \
      % (path4_5.adsorbK(), path4_5.desorbK(), path4_5.deltaG()))

path5_6 = adsorption.adsorption(s5, s6, T, P_O2, 'CO')
print("For step five:\tk(ads) = %.3e,\tk(des) = %.3e,\tDeltaG = %.3f\n" \
      % (path5_6.adsorbK(), path5_6.desorbK(), path5_6.deltaG()))

path6_7 = reactions.reactions(s6, ts6_7, s7, T)
print("For step six:\tk(forwards) = %.3e,\tk(reverse) = %.3e,\n"%(path6_7.forwardK,path6_7.reverseK))

path7_8 = adsorption.adsorption(s8, s7, T, P_CO2, 'CO2')
print("For step seven:\tk(des) = %.3e,\tk(ads) = %.3e,\tDeltaG = %.3f\n" \
      % (path7_8.desorbK(), path7_8.adsorbK(), - path7_8.deltaG()))

reactionTuple = (path1_2, path2_3, path3_4, path4_5, path5_6, path6_7, path7_8)
k_forward = (path1_2.adsorbK(), path2_3.forwardK, path3_4.desorbK(), path4_5.adsorbK(),\
    path5_6.adsorbK(), path6_7.forwardK, path7_8.desorbK())
k_reverse = (path1_2.desorbK(), path2_3.reverseK, path3_4.adsorbK(), path4_5.desorbK(),\
    path5_6.desorbK(), path6_7.reverseK, path7_8.adsorbK())
# DO lattice KMC
CeO2_stepD = lattice.lattice(latticeSize)

print(CeO2_stepD.demension)

CeO2_stepD.initialElements = speciesTuple[0].kind

print(CeO2_stepD.initialElements)

CeO2_stepD.initializeElements()


for i in range(10000):                              # MC cycle number
    k_rate = k_rate1 = 0                          # initialize k
    theReactionSite = []
    for j in range(CeO2_stepD.demension[0]):       # x direction
        for k in range(CeO2_stepD.demension[1]):   # y direction
            IS = CeO2_stepD.elements[k][j]         # find out all the IS on surface
            order = speciesTuple_kind.index(IS)    # find the IS's corresponding order
            if order == 7:                         # s8 = s1
                order = 0
            k_rate += k_forward[order] + k_reverse[order-1]  # get K_tot

    k_random = random.random()                 # random a k_random
    print(k_random,k_rate)
    for j in range(CeO2_stepD.demension[0]):       # x direction
        for k in range(CeO2_stepD.demension[1]):   # y direction
            IS = CeO2_stepD.elements[k][j]         # find out all the IS on surface
            order = speciesTuple_kind.index(IS)
            if order == 7:
                order = 0

            k_rate1 += k_forward[order]
            if  k_rate1 < k_random * k_rate:                  # once random > tot,
                pass
            elif k_rate1 > k_random * k_rate:                                  # get the reaction site
                theReactionSite.append(j)
                theReactionSite.append(k)
                theReactionSite.append(1)          # get the reaction direction
                print("break")
                break

            k_rate1 += k_reverse[order - 1]
            if  k_rate1 < k_random * k_rate:
                pass
            elif k_rate1 > k_random * k_rate:
                theReactionSite.append(j)
                theReactionSite.append(k)
                theReactionSite.append(-1)
                print("break")
                break
        else:
            continue
        break
    print(theReactionSite)
    FS_order = order + theReactionSite[2]
    if FS_order == 7:
        FS_order = 0
    FS = speciesTuple[FS_order].kind
    CeO2_stepD.elements = ((theReactionSite[0], theReactionSite[1]), IS, FS)
print(CeO2_stepD.elements)
print(CeO2_stepD.elements[0].count(s1.kind))
print(CeO2_stepD.elements[0].count(s2.kind))
print(CeO2_stepD.elements[0].count(s3.kind))
print(CeO2_stepD.elements[0].count(s4.kind))
print(CeO2_stepD.elements[0].count(s5.kind))
print(CeO2_stepD.elements[0].count(s6.kind))
print(CeO2_stepD.elements[0].count(s7.kind))
'''

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

for i in range(kmc.LoopNum):                       # MC cycle number
    k_rate = k_rate1 = 0                           # initialize k
    theReactionSite = []
    for j in range(lat.demension[0]):              # x direction
        for k in range(lat.demension[1]):          # y direction
            IS = lat.elements[k][j]                # find out all the IS on surface
            order = speciesTuple_kind.index(IS)    # find the IS's corresponding order
            if order == 7:                         # s8 = s1
                order = 0
            k_rate += k_forward[order] + k_reverse[order-1]  # get K_tot

    k_random = random.random()                 # random a k_random
    #print(k_random,k_rate)
    for j in range(lat.demension[0]):       # x direction
        for k in range(lat.demension[1]):   # y direction
            IS = lat.elements[k][j]         # find out all the IS on surface
            order = speciesTuple_kind.index(IS)
            if order == 7:
                order = 0

            k_rate1 += k_forward[order]
            if  k_rate1 < k_random * k_rate:                  # once random > tot,
                pass
            elif k_rate1 > k_random * k_rate:                                  # get the reaction site
                theReactionSite.append(j)
                theReactionSite.append(k)
                theReactionSite.append(1)          # get the reaction direction
                #print("break")
                break

            k_rate1 += k_reverse[order - 1]
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
    FS_order = order + theReactionSite[2]
    if FS_order == 7:
        FS_order = 0
    FS = speciesTuple[FS_order].kind
    lat.elements = ((theReactionSite[0], theReactionSite[1]), IS, FS)
#print(lat.elements)
for i in range(len(speciesTuple_kind)):
    print('the coverage for species {:d} = {:.3f}'.\
          format(i, float(lat.elements[0].count(speciesTuple_kind[i])) / \
            int(lat.demension[0]) * int(lat.demension[1])))


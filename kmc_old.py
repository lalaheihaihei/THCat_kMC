# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: kmc_old.py
@time: 12/1/2016 10:24 AM
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


#OLD TEST CODES for kMC, just keep here for comparation:

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

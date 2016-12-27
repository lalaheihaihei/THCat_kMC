# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: loop.py
@time: 2016/12/18 12:35
class controls the kmc loop. Need to transfer to fortran/C/C++/Cython further to accelerate loop.
"""

import random
import math
from SF import kmc
from SF import count_reaction

class Loop(object):

    def __init__(self, kmc, rate_const_dict, lat, num_of_avail_sites, time = [0], product = 0):
        """
        :param kmc: input parameters instance
        :param rate_const_dict: all reaction rate constants
        :param lat: lattice list
        :param num_of_avail_sites:  {'-6': [['-6', 3, 3]], '5': [], '-2': []...}
        :param time: list include accumulated simulation time for evergy step
        :param product: number of product
        """
        self._rate_const_dict = rate_const_dict
        self._lat = lat
        self._num_of_avail_sites = num_of_avail_sites
        self._kmc = kmc
        self._loop_n = kmc.LoopNum
        self._reactions = kmc.reactions
        self._time = time
        self._product = product
        self._coverage_list = list(map(lambda x: 0, self._kmc.kinds)) # initialize the coverage list

    def set_accum_rate(self):
        """
        generate a accum_rate dict {'-0':k_accum, '0':k_accum, '-1':k_accum, ...}
        for random select reaction in next step.

        :return: accum_rate dict {'-0':k_accum, '0':k_accum, '-1':k_accum, ...}
        """
        accum_rate = {}
        for i in self._num_of_avail_sites.keys():
            k_accum = len(self._num_of_avail_sites[i]) * self._rate_const_dict[i]
            accum_rate[i] = k_accum
        return accum_rate

    def do_random_k(self, accum_rate):
        """
        random select the elementary reaction from accum_rate {'-0':k_accum, '0':k_accum, '-1':k_accum, ...}
        For instance: sum(k0_accum:k4_accum) < random_u * k_tot < sum(k0_accum:k5_accum), so k5 is selected.

        :param accum_rate: {'-0':k_accum, '0':k_accum, '-1':k_accum, ...}
        :return: react_k: '5'
        """
        u = random.random()
        k_sum = 0
        for i in self._num_of_avail_sites.keys():
            k_sum += accum_rate[i]
            # print(k_sum, sum(accum_rate.values()))
            if k_sum > u * sum(accum_rate.values()):
                react_k = i  # find out which reaction should be process in this kmc loop
                break
        return react_k

    def do_random_site(self, react_k):
        """
        last step has selected which reaction should be done.
        in this def. find all site can do this reacion, and random select one.

        :param react_k: '5'
        :return: ['5', 3, 3], where '5' is 5th reaction. and 3, 3 means reaction site.
        """
        order_of_the_react_site = random.randint(0, len(self._num_of_avail_sites[react_k])-1)
        # print(order_of_the_react_site, self._num_of_avail_sites[react_k][order_of_the_react_site])
        return self._num_of_avail_sites[react_k][order_of_the_react_site]

    def update_time(self, accum_rate):
        """
        calculate the time go in this kmc loop.
        see Computer Physics Communications 185 (2014) 2138. Figure 2.

        :param accum_rate: the tot rate constants information
        :return: d_t, time span for this kmc step.
        """
        t_random = random.random()
        d_t = (1/sum(accum_rate.values())) * math.log(1/t_random)
        self._time.append(self._time[-1] + d_t)
        return d_t

    def update_coverage(self, d_t):
        """
        Add a time dependent coverage: see Top Catal (2014) 57:159–170 equation (6).
        self._kmc.not_count_cover could define some surface species that not involved in coverage calculation.

        :param d_t: the time span for this kmc step.
        :return: None, but update self._coverage_list
        """
        for i in range(len(self._kmc.kinds)):
            if self._kmc.kinds[i] not in self._kmc.not_count_cover:
                for j in self._lat:
                    self._coverage_list[i] += d_t * j.count(self._kmc.kinds[i])
        return None

    def update_lat(self, react_site):
        """
         self._reactions like : [[['1 ', ' CO'], ['2']], [['2 ', ' 0'], ['ts1'], ['3 ', ' 4']]
         self._reactions[int(react[1])]: select the reaction;
         self._reactions[int(react[1])][-1]: select the kinds involved in reverse reaction;
         self._reactions[int(react[1])][-1][0]: select the first kinds involved in reverse reaction;
         self._reactions[int(react[1])][-1][1]: select the second kinds involved in reverse reaction;
         self._reactions[int(react[1])][0][0]:  select the first kinds to be;
         self._reactions[int(react[1])][0][1]:  select the second kinds to be

         if the first and last element in lat is involved in reactions.(for periodic calculation.)
         the corresponding site should be update with the lat[0] or lat[-1]
         else: the lat[0] or lat[-1] should be update with the corresponding sites in real lat.

        :param react_site: ["-6", 5, 6] means the reaction "-6" was selected and it processed on site 5 and 6.
        :return: None, but update self._lat
        """
        whether_reverse = [ i for i in react_site[0]]   # split react_site[0] : "-0" to ["-","0"]
        if whether_reverse[0] == "-":
            react = react_site[0][1:]
            #print("RE@@@@@@@@@@@@@@@@", react_site, int(react), self._reactions[int(react)][0][-1][0].strip())
            if len(react_site) == 2:  # if the reaction involve one surface site
                self._lat[react_site[1][0]][react_site[1][1]] = self._reactions[int(react)][0][0][0].strip()
            if len(react_site) == 3:  # if the reaction involve two surface site
                self._lat[react_site[1][0]][react_site[1][1]] = self._reactions[int(react)][0][0][0].strip()
                self._lat[react_site[2][0]][react_site[2][1]] = self._reactions[int(react)][0][0][1].strip()
            if len(react_site) == 4:  # if the reaction involve three surface site
                self._lat[react_site[1][0]][react_site[1][1]] = self._reactions[int(react)][0][0][0].strip()
                self._lat[react_site[2][0]][react_site[2][1]] = self._reactions[int(react)][0][0][1].strip()
                self._lat[react_site[3][0]][react_site[3][1]] = self._reactions[int(react)][0][0][2].strip()
        else:
            react = react_site[0]
            #print("@@@@@@@@@@@@@@@@", react_site, int(react), self._reactions[int(react)][0][-1][0].strip())
            if len(react_site) == 2:  # if the reaction involve one surface site
                self._lat[react_site[1][0]][react_site[1][1]] = self._reactions[int(react)][0][-1][0].strip()
            if len(react_site) == 3:  # if the reaction involve two surface site
                self._lat[react_site[1][0]][react_site[1][1]] = self._reactions[int(react)][0][-1][0].strip()
                self._lat[react_site[2][0]][react_site[2][1]] = self._reactions[int(react)][0][-1][1].strip()
            if len(react_site) == 4:  # if the reaction involve three surface site
                self._lat[react_site[1][0]][react_site[1][1]] = self._reactions[int(react)][0][-1][0].strip()
                self._lat[react_site[2][0]][react_site[2][1]] = self._reactions[int(react)][0][-1][1].strip()
                self._lat[react_site[3][0]][react_site[3][1]] = self._reactions[int(react)][0][-1][2].strip()
        '''
        if self._kmc.periodic == "1":
            if 0 in react_site:
                self._lat[-2] = self._lat[0]
                self._lat[-1] = self._lat[1]
            elif len(self._lat)-1 in react_site:
                self._lat[0] = self._lat[-2]
                self._lat[1] = self._lat[-1]
            else:
                self._lat[0] = self._lat[-2]
                self._lat[-1] = self._lat[1]
        '''
        return None

    def update_num_of_avail_sites(self):
        """
        this step is quite slow because of re-evaluate the n_avail by search whole lattice.
        need to be optimized for large lattice, see: Computer Physics Communications 185 (2014) 2138.
        :return: None, but update self._num_of_avail_sites
        """
        self._num_of_avail_sites = count_reaction.initialize_num_of_avail_sites(self._kmc, self._lat)
        return None

    def count_product(self, react_site):
        """
        This def counts product number. if the selected reaction is same as input files count key (kmc.count_product)
        count +1, on the contrary, if the reverse reaction is done. count -1

        :param react_site: ["-6", 5, 6] means the reaction "-6" was selected and it processed on site 5 and 6.
        :return: None, count product number.
        """
        for i in self._kmc.count_product:
            if react_site[0] == i:
                self._product += 1
                #print("+1")
            elif react_site[0] == "-"+i:
                self._product -= 1
                #print("-1")

    def do_kmc_loop(self):
        print("\n&&&&&&&&&&&&&&&&&&&&&&&&&   prepare kMC loop   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n")
        #print("rate_const_dict :", self._rate_const_dict)
        #print("The lattice :", self._lat)
        for (k, v) in self._num_of_avail_sites.items():
            if v != []:
                print("reaction %s:" % k, v)
        # print("number of the available sites for specific elementary reaction is: \n", self._num_of_avail_sites)
        print("loop number is", self._loop_n)
        print("\n&&&&&&&&&&&&&&&&&&&&&&&&&      start loop      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n")

        for i in range(self._loop_n):

            # update accumulate rate for every elementary reaction step.
            accum_rate = self.set_accum_rate()

            # update time
            d_t = self.update_time(accum_rate)

            # store coverage infomation
            self.update_coverage(d_t)

            # select the reaction
            react_k = self.do_random_k(accum_rate)

            # select the lattice site
            react_site = self.do_random_site(react_k)

            # update surface species
            self.update_lat(react_site)

            # update available reactions dictionary
            self.update_num_of_avail_sites()

            # count product number
            self.count_product(react_site)

            # print for every 1000 steps
            if i%1000 == 0:
                print("\n################  loop %d  ####################\n" % i)
                for (k, v) in accum_rate.items():
                    if v != 0:
                        print("accumulated rate %3s : %.3e" % (k, v))
                print("total time for loop 0 -%8d is %.3e:" % (i, self._time[-1]))
                print("kMC step %s, reaction %s happens on lattice site:" % (i, react_k), react_site)
                print("lattice model:")
                for i in self._lat:
                    print(i)
                print("product number is %d", self._product)


        '''
        print("number of SAs:", int(self._kmc.num_SA)-1)
        print("there are %d products in %.2f s" % (self._product, self._time[-1]))
        # self.kmc.num_SA-1 is the number of SAs because of periodic condition cutoff the first lattice site.
        # print("the real active sites number is %d - %d = %d" % \
        #      (int(self._kmc.num_SA), self._count_cut_num_of_active, int(self._kmc.num_SA) - self._count_cut_num_of_active))
        # print("tof = ", math.log(self._product/ (self._time[-1] * (int(self._kmc.num_SA)-self._count_cut_num_of_active)) ))

        '''
        #self._coverage_list = list(map(lambda x: x/sum(self._coverage_list), self._coverage_list))
        #for i in range(len(self._coverage_list)):
        #    print("coverage of surface species %s is %.3f:" % (self._kmc.kinds[i], self._coverage_list[i]))

        return None

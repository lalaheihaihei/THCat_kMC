# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: loop.py
@time: 2016/12/18 12:35
"""

import random
import math
from oneD_SAmodule import kmc

class Loop(object):

    def __init__(self, kmc, rate_const_dict, lat, num_of_avail_sites, count_cut_num_of_active, time = [0], product = 0):
        self._rate_const_dict = rate_const_dict
        self._lat = lat
        self._num_of_avail_sites = num_of_avail_sites
        self._kmc = kmc
        self._loop_n = kmc.LoopNum
        self._reactions = kmc.reactions
        self._count_cut_num_of_active = count_cut_num_of_active
        self._time = time
        self._product = product
        self._coverage_list = list(map(lambda x: 0, self._kmc.kinds)) # initialize the coverage list

    def set_accum_rate(self):
        accum_rate = {}
        for i in self._num_of_avail_sites.keys():
            k_accum = len(self._num_of_avail_sites[i]) * self._rate_const_dict[i]
            accum_rate[i] = k_accum
        return accum_rate

    def do_random_k(self, accum_rate):
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
        order_of_the_react_site = random.randint(0, len(self._num_of_avail_sites[react_k])-1)
        # print(order_of_the_react_site, self._num_of_avail_sites[react_k][order_of_the_react_site])
        return self._num_of_avail_sites[react_k][order_of_the_react_site]

    def update_time(self, accum_rate):
        t_random = random.random()
        d_t = (1/sum(accum_rate.values())) * math.log(1/t_random)
        self._time.append(self._time[-1] + d_t)
        return d_t

    def update_coverage(self, d_t):
        # add a time dependent coverage: see Top Catal (2014) 57:159–170 equation (6).
        for i in range(len(self._kmc.kinds)):
            if self._kmc.kinds[i] not in self._kmc.not_count_cover:
                self._coverage_list[i] += d_t * self._lat.count(self._kmc.kinds[i])
        return None

    def update_lat(self, react_site):
        react = [ i for i in react_site[0]]   # split react_site[0] : "-0" to ["-","0"]
        if react[0] == "-":
            """self._reactions like : [[['1 ', ' CO'], ['2']], [['2 ', ' 0'], ['ts1'], ['3 ', ' 4']]
             self._reactions[int(react[1])]: select the reaction;
             self._reactions[int(react[1])][-1]: select the kinds involved in reverse reaction;
             self._reactions[int(react[1])][-1][0]: select the first kinds involved in reverse reaction;
             self._reactions[int(react[1])][-1][1]: select the second kinds involved in reverse reaction;
             self._reactions[int(react[1])][0][0]:  select the first kinds to be;
             self._reactions[int(react[1])][0][1]:  select the second kinds to be"""
            self._lat[react_site[1]] = self._reactions[int(react[1])][0][0].strip()
            if react_site[1] != react_site[2]:  # if the reaction involve two surface site
                self._lat[react_site[2]] = self._reactions[int(react[1])][0][1].strip()
        else:
            self._lat[react_site[1]] = self._reactions[int(react[0])][-1][0].strip()
            if react_site[1] != react_site[2]:  # if the reaction involve two surface site
                self._lat[react_site[2]] = self._reactions[int(react[0])][-1][1].strip()
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
        return None

    def update_num_of_avail_sites(self):
        """
        this step is quite slow because of re-evaluate the n_avail by search whole lattice.
        :return:
        """
        self._num_of_avail_sites = kmc.initialize_num_of_avail_sites(self._kmc, self._lat)
        return None

    def count_product(self, react_site):
        if react_site[0] == self._kmc.count_product:
            self._product += 1
            #print("+1")
        elif react_site[0] == "-"+self._kmc.count_product:
            self._product -= 1
            #print("-1")


    def do_kmc_loop(self):
        print("\n\n&&&&&&&&&&&&&&&&&&&&&&&&&   prepare kMC loop   &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        print("rate_const_dict :", self._rate_const_dict)
        print("The lattice :", self._lat)
        print("number of the available sites for specific elementary reaction is: \n", self._num_of_avail_sites)
        print("loop number is", self._loop_n)

        print("&&&&&&&&&&&&&&&&&&&&&&&&&      start loop      &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")

        for i in range(self._loop_n):
            accum_rate = self.set_accum_rate()
            # print("accumulated rate is", accum_rate)
            d_t = self.update_time(accum_rate)
            self.update_coverage(d_t)
            # print(self._time)

            react_k = self.do_random_k(accum_rate)

            react_site = self.do_random_site(react_k)
            # print("kMC step %d, reaction %s is processed on lattice site %d %d." % (i, react_k, react_site[1], react_site[2]))
            self.update_lat(react_site)
            # print(react_site, self._lat)
            self.update_num_of_avail_sites()
            # print(self._num_of_avail_sites, "\n")
            self.count_product(react_site)
        print("lattice:", self._lat)
        print("number of SAs:", int(self._kmc.num_SA)-1)
        print("there are %d products in %.2f s" % (self._product, self._time[-1]))
        # self.kmc.num_SA-1 is the number of SAs because of periodic condition cutoff the first lattice site.
        print("the real active sites number is %d - %d = %d" % \
              (int(self._kmc.num_SA), self._count_cut_num_of_active, int(self._kmc.num_SA) - self._count_cut_num_of_active))
        print("tof = ", math.log(self._product/ (self._time[-1] * (int(self._kmc.num_SA)-self._count_cut_num_of_active)) ))
        self._coverage_list = list(map(lambda x: x/sum(self._coverage_list), self._coverage_list))
        for i in range(len(self._coverage_list)):
            print("coverage of surface species %s is %.3f:" % (self._kmc.kinds[i], self._coverage_list[i]))
        return None
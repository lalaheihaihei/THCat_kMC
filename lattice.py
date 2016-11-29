# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: lattice.py
@time: 11/29/2016 11:27 AM
A lattice class control the initialization and evolution of surface species.
"""

#import species

class lattice(object):


    def __init__(self,demension):
        self._demension = demension


    # get the demension on surface lattice: demension[0]*demension[1]
    @property
    def demension(self):
        return self._demension

    @demension.setter
    def demension(self,d):
        if self.test2Dtuple(d):
            self._demension = d


    # initialize the adsorbed species on surface
    @property
    def initialElements(self):
        return self._initialElements

    @initialElements.setter
    def initialElements(self,a):
        self._initialElements = [[i for i in range(self._demension[0])]\
                                 for i in range(self._demension[1])]
        for i in range(0, self._demension[0]):
            for j in range(0, self._demension[1]):
                self._initialElements[j][i] = a
        return self._initialElements


    # initialize self._elements
    def initializeElements(self):
        self._elements = self._initialElements
        pass


    # evolution of self._elements
    @property
    def elements(self):
        return self._elements

    @elements.setter
    def elements(self, siteAndISFS):
        site = siteAndISFS[0]
        IS = siteAndISFS[1]
        FS = siteAndISFS[2]
        if self.test2Dtuple(site):
            pass
        if not isinstance(IS,str) and isinstance(FS,str):
            raise ValueError('IS and FS must be a spieces of str')
        if self._elements[site[1]][site[0]] != IS:
            raise IOError('Error: the Initial state is not same as the one in dict')
        self._elements[site[1]][site[0]] = FS
        return self._elements




    # Test all the input 2D tuples in this class.
    def test2Dtuple(self,d):
        if not isinstance(d,tuple):
            raise ValueError('demension should be a tuple')
        elif len(d) != 2:
            raise ValueError('demension should be a tuple with 2 elements')
        elif not isinstance(d[0],int) and isinstance(d[1],int):
            raise ValueError('demension\'s elements should be int')
        else:
            return True



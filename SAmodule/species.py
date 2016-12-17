# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: species.py
@time: 11/29/2016 10:52 AM
A species class, including TS, which all contain frequency and energy information.
"""


class species(object):

    def __init__(self, kind = 'CO', freq = (100,1000,2000), energy = 10):
        self._kind = kind
        self._freq = freq
        self._energy = energy

    # set the adsorbed speices.
    @property
    def kind(self):
        return self._kind

    @kind.setter
    def kind(self,a):
        if not isinstance(a,str):
            raise ValueError('kind must be string')
        self._kind = a

    # set the frequency.
    @property
    def freq(self):
        return self._freq

    @freq.setter
    def freq(self,values):
        if not isinstance(values,tuple):
            raise ValueError('freq must be tuple')
        for i in values:
            if not isinstance(i,float):
                raise ValueError('freq tuple\' elements must be float')
        self._freq = values

    # set the energy.
    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self,value):
        if not isinstance(value,float):
            raise isinstance('energy must be float')
        self._energy = value

    # get the ZPE energy.
    @property
    def ZPE(self):
        return sum(list(map(lambda x: 1/2 * x * 1.23981e-4,self._freq)))

    @property
    def correctEnergy(self):
        # print('ZPEZPEZPEZPEZPE', self._energy - self.ZPE) # for bug fix.
        return self._energy + self.ZPE   # BUG FIXED: before Dec 18th, the ZPE is minus...



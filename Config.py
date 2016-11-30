# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: Config.py.py
@time: 11/30/2016 11:38 AM
Read configuration of the species, reactions, lattice and so on
from a config file. This function is a little stupid, but it works.
"""

import configparser, os

class Parameters(object):
    """
    A class to read and write the config file in THCat-kMC.
    Using a class to handle a bunch of parameters might be a good idea.
    """


    def __init__(self, FileName = 'config.txt'):
        """
        Default value for parameters
        """
        # Set some necessary directories.
        self.HomeDir = os.getcwd()
        config = configparser.ConfigParser()
        config.read(FileName)
        print(config.sections())

        # Set some common parameters for kMC.
        secCommon = config['common']
        self.RunType = secCommon['RunType']
        print('Reading configuration file for %s calculation' % (self.RunType))

        self.Temperature = secCommon['Temperature']
        print('kMC will be performed at %s K' % (self.Temperature))

        self.gas = secCommon['Gas'].split(' ')
        print('The gas molecules involved in the system including:', self.gas)

        self.GasPressure = secCommon['GasPressure'].split(' ')

        # Set some common parameters for kMC.
        secSAspecies = config['SAspecies']
        self.kinds = tuple(map(lambda x : x.strip(), secSAspecies['kinds'].split(',')))
        print(self.kinds)

        self.ifKindsTS = tuple(map(lambda x : int(x), secSAspecies['ifKindsTS'].split(',')))
        print(self.ifKindsTS)

        self.kindsFreq = list(map(lambda x : x.strip().split(' '), secSAspecies['kindsFreq'].split(',')))
        for i in range(len(self.kindsFreq)):
            self.kindsFreq[i] = list(map(lambda x : float(x), self.kindsFreq[i]))
        print(self.kindsFreq)


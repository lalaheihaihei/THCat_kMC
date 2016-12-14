# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: Config.py.py
@time: 11/30/2016 11:38 AM
Read configuration of the species, reactions, lattice and so on
from a config file. This function is a little stupid, but it works.
"""

import configparser
import os
import time


class Parameters(object):
    """
    A class to read and write the config file in THCat-kMC.
    Using a class to handle a bunch of parameters might be a good idea.
    """

    def __init__(self, filename='config-N2.txt'):
        """
        Default value for parameters
        """
        # Set necessary directories and job time.
        self.__filename = filename
        self.HomeDir = os.getcwd()
        Configfilepath = ('/'.join([self.HomeDir, filename]))
        config = configparser.ConfigParser()
        config.read(Configfilepath)
        starttime = time.clock()
        print('Starting calculation at')
        print(time.strftime("%H:%M:%S on %a %d %b %Y"))


        # Set some common parameters for kMC.
        secCommon = config['common']

        self.Title = secCommon['Title']
        print('THCat_kMC Start! For system: %s ' % (self.Title))
        print('\'Config.txt\' include %d parts: ' % (len(config.sections())), config.sections(), '\n')

        self.RunType = secCommon['RunType']
        print('Reading configuration file for %s calculation' % (self.RunType))

        self.Temperature = float(secCommon['Temperature'])
        print('kMC will be performed at %s K' % (self.Temperature))

        self.gas = secCommon['Gas'].split(' ')
        print('The gas molecules involved in the system including:', self.gas)

        self.GasPressure = tuple(map(lambda x : float(x), secCommon['GasPressure'].split(',')))

        self.latticeSize = tuple(map(lambda x : int(x), secCommon['latticeSize'].split(',')))

        self.LoopNum = int(secCommon['LoopNum'].strip())
        print('Loop number of kMC: %d' % (self.LoopNum))



        # Set some SA_kMC parameters for kMC.
        if self.RunType == 'SA':
            sacSA = config['SA']

            self.kinds = tuple(map(lambda x : x.strip(), sacSA['kinds'].split(',')))
            print(self.kinds)

            self.isKindsTS = tuple(map(lambda x : int(x), sacSA['isKindsTS'].split(',')))
            print(self.isKindsTS)

            self.kindsEnergy = tuple(map(lambda x : float(x.strip()), sacSA['kindsEnergy'].split(',')))
            print(self.kindsEnergy)

            self.kindsFreq = list(map(lambda x: x.strip().split(' '), sacSA['kindsFreq'].split(',')))
            for i in range(len(self.kindsFreq)):
                self.kindsFreq[i] = list(map(lambda x: float(x), self.kindsFreq[i]))
            print(self.kindsFreq)

            self.reactions = list(map(lambda x: x.strip().split('<-->'), sacSA['reactions'].split(',')))
            for i in range(len(self.reactions)):
                self.reactions[i] = list(map(lambda x: x.strip().split('+'), self.reactions[i]))
            print(self.reactions)

            self.reactionsKind = tuple(map(lambda x : str(x.strip()), sacSA['reactionsKind'].split(',')))
            print(self.reactionsKind)


        # Set some Surface_kMC parameters for kMC.
        if self.RunType == 'SURFACE':
            sacSF = config['SF']
            self.kinds = tuple(map(lambda x : x.strip(), sacSF['kinds'].split(',')))
            print(self.kinds)

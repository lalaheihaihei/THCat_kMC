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

    def __init__(self, filename='config-TiO2.txt'):
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

        # Set some 1D_SA_kMC parameters for kMC.
        if self.RunType == '1D_SA':
            print("###########################    1D_SA module    ########################################")
            sac1D_SA = config['1D_SA']

            self.num_SA = sac1D_SA['number_SA'].strip()

            self.kinds = tuple(map(lambda x : x.strip(), sac1D_SA['kinds'].split(',')))
            print('surface species contain:', self.kinds)
            self.TSkinds = tuple(map(lambda x: x.strip(), sac1D_SA['TSkinds'].split(',')))
            print("surface TS species contain:", self.TSkinds)

            self.reactions = list(map(lambda x: x.strip().split('<-->'), sac1D_SA['reactions'].split(',')))
            for i in range(len(self.reactions)):
                self.reactions[i] = list(map(lambda x: x.strip().split('+'), self.reactions[i]))
            print("All elementary steps contain:", self.reactions)

            self.Energies = list(map(lambda x : x.strip().split(',') , sac1D_SA['Energies'].split(';')))
            for i in range(len(self.Energies)):
                self.Energies[i] = list(map(lambda x: float(x.strip()), self.Energies[i]))
            self.Energies = tuple(self.Energies)
            print("energy for all elementary steps" ,self.Energies)

            self.Freq = list(map(lambda x: x.strip().split(' '), sac1D_SA['Freq'].split(',')))
            for i in range(len(self.Freq)):
                self.Freq[i] = list(map(lambda x: float(x), self.Freq[i]))
            print("All elementary steps' freq is ,", self.Freq)

            self.reactionsKind = tuple(map(lambda x : str(x.strip()), sac1D_SA['reactionsKind'].split(',')))
            print("Kinds of elementary steps: for calculation of reaction rate K:", self.reactionsKind)

            self.count_product = sac1D_SA['countProduct'].strip()

            self.not_count_cover = tuple(map(lambda x : str(x.strip()), sac1D_SA['notCountCover'].split(',')))

            self.periodic = sac1D_SA['periodic'].strip()

        # Set some Surface_kMC parameters for kMC.
        if self.RunType == 'SF':
            print("###########################    SF module    ########################################")
            SF = config['SF']

            self.num_SA = SF['number_SA'].strip()

            self.kinds = SF['kinds'].split(',')
            if len(self.kinds) == 1:
                self.kinds = self.kinds[0].split('-')
                self.kinds = tuple(map(lambda x: int(x.strip()), self.kinds))
                self.kinds = [str(i) for i in range(self.kinds[0], self.kinds[1]+1)]
                print('surface species contain:', self.kinds)
            else:
                self.kinds = tuple(map(lambda x: x.strip(), SF['kinds'].split(',')))
                print('surface species contain:', self.kinds)

            self.reactions = list(map(lambda x: x.strip(), SF['reactions'].split('\n')))
            for i in range(len(self.reactions)):
                self.reactions[i] = self.reactions[i].split(";")
                self.reactions[i][0] = self.reactions[i][0].split('<-->')
                for j in range(len(self.reactions[i][0])):
                    self.reactions[i][0][j] = list(map(lambda x: x.strip(), self.reactions[i][0][j].split('+')))
                self.reactions[i][1] = self.reactions[i][1].strip()
                self.reactions[i][2] = list(map(lambda x: float(x.strip()), self.reactions[i][2].split(",")))


            #    self.reactions[i] = list(map(lambda x: x.split(';'), self.reactions[i]))
            print("All elementary steps contain:", self.reactions)


'''
            self.Energies = list(map(lambda x: x.strip().split(','), SF['Energies'].split(';')))
            for i in range(len(self.Energies)):
                self.Energies[i] = list(map(lambda x: float(x.strip()), self.Energies[i]))
            self.Energies = tuple(self.Energies)
            print("energy for all elementary steps", self.Energies)

            self.Freq = list(map(lambda x: x.strip().split(' '), SF['Freq'].split(',')))
            for i in range(len(self.Freq)):
                self.Freq[i] = list(map(lambda x: float(x), self.Freq[i]))
            print("All elementary steps' freq is ,", self.Freq)

            self.reactionsKind = tuple(map(lambda x: str(x.strip()), SF['reactionsKind'].split(',')))
            print("Kinds of elementary steps: for calculation of reaction rate K:", self.reactionsKind)

            self.count_product = SF['countProduct'].strip()

            self.not_count_cover = tuple(map(lambda x: str(x.strip()), SF['notCountCover'].split(',')))

            self.periodic = SF['periodic'].strip()
'''

def not_empty(self, s):
    return s and s.strip()
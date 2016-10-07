# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 20:19:49 2016

@author: aschmitt
"""


import os
import sys
import re
import random

import matplotlib as mpl
mpl.use('Qt5Agg')

import matplotlib.pyplot as plt
import numpy as np

import argparse

class Reference(object):
    def __init__(self, referencepath, rtype="pfs"):
        self.logTagStr = "Reference"
        
        self.name = os.path.basename(referencepath)
        self.rtype = rtype
        self.referencepath = referencepath
        self.n = 0

        self.ids = []
        self.redshifts = []
        self.mags = []
        self.flags = []
        self.sfrs = []
        self.ebmvs = []
        self.sigmas = []
        self.loghalphas= []

        # idx variables, default for VVDS
        self.nRefValues = 7
        self.iRefID = 0
        self.iRefZ = 4
        self.iRefMag = 6
        self.iRefFlag = 5
        self.iRefSFR = -1
        self.iRefEBmV = -1
        self.iRefSigma = -1
        self.iRefLogHalpha = -1        
        
        if self.rtype=="pfs":
            self.setPFSRefFileType()
            
        self.load()

    def setMuseRefFileType(self):
        self.nRefValues = 2
        self.iRefZ = 1
        self.iRefMag = -1
        self.iRefFlag = -1
        
        self.iRefSFR = -1
        self.iRefEBmV = -1
        self.iRefSigma = -1
        
    def setVUDSRefFileType(self):
        self.nRefValues = 6
        self.iRefZ = 3
        self.iRefMag = 5
        self.iRefFlag = 4  
        
    def setVVDSRefFileType(self):
        self.nRefValues = 7
        self.iRefZ = 4
        self.iRefMag = 6
        self.iRefFlag = 5
    
    def setVVDS2RefFileType(self):
        self.nRefValues = 7
        self.iRefZ = 5
        self.iRefMag = 4
        self.iRefFlag = 6
    
    def setPFSRefFileType(self):
        self.nRefValues = 7
        self.iRefZ = 1
        self.iRefMag = 2
        self.iRefFlag = -1
        self.iRefSFR = 5
        self.iRefEBmV = 4
        self.iRefSigma = 6
        
    def setKeckRefFileType(self):
        """
        This corresponds to the Keck scientists reference files, ex. refz_vlebrun_cos-m4n1.txt
        """
        self.nRefValues = 8
        self.iRefID = 3
        self.iRefZ = 1
        self.iRefMag = 2
        self.iRefFlag = 7
        self.iRefSFR = 5
        self.iRefEBmV = 4
        self.iRefSigma = 6
    
    def setSIMULMRefFileType(self):
        self.nRefValues = 8
        self.iRefZ = 1
        self.iRefMag = 2
        self.iRefFlag = -1
        self.iRefSFR = 5
        self.iRefEBmV = 4
        self.iRefISM = 6
        self.iRefSigma = 7
        
    def setSIMUEuclid2016RefFileType(self):
        self.nRefValues = 8
        self.iRefZ = 1
        self.iRefMag = 2
        self.iRefFlag = -1
        self.iRefSFR = 5
        self.iRefEBmV = 4
        self.iRefISM = -1
        self.iRefSigma = 6
        self.iRefLogHalpha = 7    
        
    def setAriadneRefFileType(self):
        self.nRefValues = 2
        self.iRefZ = 1
        self.iRefMag = -1
        self.iRefFlag = -1
        self.iRefSFR = -1
        self.iRefEBmV = -1
        self.iRefISM = -1
        self.iRefSigma = -1
        self.iRefLogHalpha = -1   
        
    def setNoRefFileType(self):
        self.nRefValues = 0
        self.iRefZ = -1
        self.iRefMag = -1
        self.iRefFlag = -1
        self.iRefSFR = -1
        self.iRefEBmV = -1
        self.iRefISM = -1
        self.iRefSigma = -1
        self.iRefLogHalpha = -1
            

    def load(self):
        """
        load the ref file data
        """ 
        print("loadref with type={}".format(self.rtype))

        f = open(self.referencepath)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                data = lineStr.split("\t")
                if len(data)<2:
                    data = lineStr.split(" ")
                data = [r for r in data if r != '']
                
                _id = -1
                if self.iRefID!=-1:
                    _id = str(data[self.iRefID])
                self.ids.append(_id)
                
                _redshift = -1
                if self.iRefZ!=-1:
                    _redshift = float(data[self.iRefZ])
                self.redshifts.append(_redshift)
                
                _mag = -1
                if self.iRefMag!=-1:
                    _mag = float(data[self.iRefMag])
                self.mags.append(_mag)
                
                _flag = -1
                if self.iRefFlag!=-1:
                    _flag = str(data[self.iRefFlag])
                self.flags.append(_flag)
                
                _sfr = -1
                if self.iRefSFR!=-1:
                    _sfr = str(data[self.iRefSFR])
                self.sfrs.append(_sfr)
                
                _ebmv = -1
                if self.iRefEBmV!=-1:
                    _ebmv = str(data[self.iRefEBmV])
                self.ebmvs.append(_ebmv)
                
                _sigma = -1
                if self.iRefSigma!=-1:
                    _sigma = str(data[self.iRefSigma])
                self.sigmas.append(_sigma)
                
                _loghalpha = -1
                if self.iRefLogHalpha!=-1:
                    _loghalpha = str(data[self.iRefLogHalpha])
                self.loghalphas.append(_loghalpha)
              
    def getSFR(self, idList):
        sfrList = []
        for _idTarget in idList:
            found=False
            for k, _idThis in enumerate(self.ids):
                if _idThis in _idTarget:
                    sfrList.append(self.sfrs[k])
                    found = True
            if found == False:
                print("COuld not find id = {}".format(_idTarget))
                stop
        return sfrList
    
    def plot(self):
        plt.figure()
        plt.plot(self.redshifts, 'x')
        plt.grid()
        plt.xlabel("#")
        plt.ylabel("reference {}".format("redshift"))
        plt.show()
        
        
def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--referencePath", dest="referencePath", default="",
                    help="path to the reference file to be loaded")    
                    
    parser.add_argument("-t", "--type", dest="type", default="pfs",
                    help="reference file type, choose between 'vvds1', 'vvds2' or 'pfs', 'vuds', 'simulm', 'simueuclid2016', 'keck', 'no'")
          
    options = parser.parse_args()
    #print(options)

    if os.path.exists(options.referencePath) :
        rpath = os.path.abspath(options.referencePath)
        print('using full path: {0}'.format(rpath))
        ref = Reference(rpath, options.type)

        
        ref.plot()
    else :
        print("Error: invalid arguments")
        exit()
    
    
    
def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print("Reference")
    Main( sys.argv )
    

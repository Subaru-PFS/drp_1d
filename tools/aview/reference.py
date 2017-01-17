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
if "DISPLAY" in os.environ:
    mpl.use('Qt5Agg')
else:
    mpl.use('Agg')
    print("WARNING: No display available, Agg backend used as fallback")

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
        self.loghalphas = []

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
        
        print("INFO: loading ref. with rtype = {}".format(self.rtype))
        if self.rtype=="pfs":
            self.setPFSRefFileType()        
        elif self.rtype=="simueuclid2016":
            self.setSIMUEuclid2016RefFileType()
        elif self.rtype=="vvds":
            self.setVVDSRefFileType()
        else:
            print("Type not understood. Please check the ref. type arg., aborting...")
            exit()
            
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
        self.iRefID = 0
        self.iRefZ = 4
        self.iRefMag = 6
        self.iRefFlag = 5
        self.iRefSFR = -1
        self.iRefEBmV = -1
        self.iRefSigma = -1
        self.iRefLogHalpha = -1 
    
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
                    _sfr = float(data[self.iRefSFR])
                self.sfrs.append(_sfr)
                
                _ebmv = -1
                if self.iRefEBmV!=-1:
                    _ebmv = float(data[self.iRefEBmV])
                self.ebmvs.append(_ebmv)
                
                _sigma = -1
                if self.iRefSigma!=-1:
                    _sigma = float(data[self.iRefSigma])
                self.sigmas.append(_sigma)
                
                _loghalpha = -1
                if self.iRefLogHalpha!=-1:
                    _loghalpha = float(data[self.iRefLogHalpha])
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
        
    def filterIdList(self, idList, magmin=-1, magmax=50, zmin=-1, zmax=50, sfrmin=-1, sfrmax=1e6):
        """
        return a list of indexes with the mag, redshift and sfr matching the filter range
        """
        print("about ot filter with zmin={}, zmax={}, magmin={}, magmax={}, sfrmin={}, sfrmax={}".format(zmin, zmax, magmin, magmax, sfrmin, sfrmax))
        indexes = []
        for k, _idTarget in enumerate(idList):
            #print("idTarget={}".format(_idTarget))
            found=False
            for kThis, _idThis in enumerate(self.ids):
                #print("_idThis={}".format(_idThis))
                if _idThis in _idTarget:
                    #print("_idThis={} matching found".format(_idThis))
                    #print("redshift={}".format(self.redshifts[kThis]))
                    #print("mag={}".format(self.mags[kThis]))
                    #print("sfr={}".format(self.sfrs[kThis]))
                    
                    if self.mags[kThis]>=magmin and self.mags[kThis]<=magmax:
                        if self.redshifts[kThis]>=zmin and self.redshifts[kThis]<=zmax:
                            #print("sfr k = {}".format(self.sfrs[kThis]))
                            #print("sfr k = {}".format(type(self.sfrs[kThis])))
                            #print("sfrmin = {}".format(sfrmin))
                            #print("sfrmax = {}".format(sfrmax))
                            #print("sfrmax = {}".format(type(sfrmax)))
                            if self.sfrs[kThis]<=sfrmax and self.sfrs[kThis]>=sfrmin:
                                indexes.append(k)
                                found = True
                    
                if found:
                    break
        
        return indexes
        
    def export(self, outputFilePath):
        """
        export in the euclidsim2016 format by default, coz it's the most exhaustive format
        """
        f = open(outputFilePath, 'w')
        f.write("#id   z_true   ref_mag_r01   spcid   ebv   sfr   sigma   logf_halpha_ext   logf_hbeta_ext   logf_o2_ext   logf_o3_ext   logf_n2_ext   logf_s2_ext")
        f.write("\n")        
        for k, _idThis in enumerate(self.ids):            
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            self.ids[k], 
            self.redshifts[k], 
            self.mags[k],
            -1,
            self.ebmvs[-1],
            self.sfrs[k], 
            self.sigmas[k],
            self.loghalphas[k],
            -1,
            -1,
            -1,
            -1,
            -1))
            
        f.close()
        
    def convertRedshift(self, initialLbda=6562.8, finalLbda=3200.0):
        for k, _id in enumerate(self.ids):            
            newRedshift = (1+self.redshifts[k])*initialLbda/finalLbda-1
            self.redshifts[k] = newRedshift
    
    def getRedshift(self, idTag):
        """
        return the redshift value for the input id (=nameTag)
        """
        index = self.findIdx(idTag)
        if index == -1:
            stop
        z = self.redshifts[index]
        return z
        
    def findIdx(self, idNameTag):
        """
        returns the idx of the ref entry corresponding to the input id/nameTag
        """
        p = [i for i,x in enumerate(self.ids) if str(x) in idNameTag] 
        if len(p) == 1:
            #print "OK : index found : ref={0} and calc={1}".format(s,p)
            return p[0]
        elif len(p) > 1:
            print "ERROR : multiple index found for a given nametag!!".format()
            print "ERROR : idNameTag = {}\n".format(idNameTag)
            return -1
        else:
            print "ERROR : index not found for idNameTag={} in the reference ids list!!".format(idNameTag)
            return -1
        
    def getTag(self, idx):
        """
        returns a string tag with obj description for plotting
        """
        str_ = "z={:.5f} mag={:.2f}, sfr={:.3f}".format(self.redshifts[idx], self.mags[idx], self.sfrs[idx])
        return str_
    
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

        enableConvertAndExport = False
        if enableConvertAndExport:
            ref.convertRedshift()
            exportFileName = "{}_convertedzforoii.txt".format(os.path.splitext(os.path.basename(rpath))[0])
            exportFilePath = os.path.join(os.path.split(rpath)[0], exportFileName)
            print("export path is {}".format(exportFilePath))
            if os.path.exists(exportFilePath):
                print("Export path already exists, aborting export...")
            else:
                ref.export(exportFilePath)
        
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
    

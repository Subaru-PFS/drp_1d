# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 20:19:49 2016

@author: aschmitt
"""
## python 3 compatibility imports ###
#from __future__ import division, print_function
#from future_builtins import *
#from builtins import range
#from future.builtins.disabled import apply, cmp, file, raw_input, xrange
## ##### # ############# ####### ###

import os
import sys
import re
import random
import inspect

import matplotlib as mpl
if "DISPLAY" in os.environ:
    #mpl.use('Qt5Agg')
    #mpl.use('Gtk3Agg')
    pass
else:
    mpl.use('Agg')
    print("WARNING: No display available, Agg backend used as fallback")

import matplotlib.pyplot as plt
import numpy as np

try:
    subfolder = "../stats"
    cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],subfolder)))
    if cmd_subfolder not in sys.path:
        #print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
        sys.path.insert(0, cmd_subfolder)
     
    import lstats
except:
    print("Import ERROR: unable to load the lstats package. some functionnalities will be unavailable...")

import argparse

class Reference(object):
    def __init__(self, referencepath, rtype="pfs", norefIdList=[]):
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
        self.avstars = []
        self.avgas = []

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
        self.iRefAvStars = -1
        self.iRefAvGas = -1        
        
        print("INFO: loading ref. with rtype = {}".format(self.rtype))
        if self.rtype=="pfs":
            self.setPFSRefFileType()        
        elif self.rtype=="simueuclid2016":
            self.setSIMUEuclid2016RefFileType()
        elif self.rtype=="simueuclid2017":
            self.setSIMUEuclid2017RefFileType()
        elif self.rtype=="pfs201703":
            self.setPFS201703RefFileType()
        elif self.rtype=="vvds":
            self.setVVDSRefFileType()
        elif self.rtype=="noref":
            self.setNoRefFileType()
        elif self.rtype=="ariadne":
            self.setAriadneRefFileType()
        elif self.rtype=="muse2":
            self.setMuse2RefFileType()
        elif self.rtype=="simple":
            self.setSimpleRefFileType()
        else:
            print("Type not understood. Please check the ref. type arg., aborting...")
            exit()
            
        if not self.rtype=="noref":
            self.load()
        else:
            self.fakeLoad(norefIdList)

    def setSimpleRefFileType(self):
        self.nRefValues = 2
        self.iRefZ = 1
        self.iRefMag = -1
        self.iRefFlag = -1
        
        self.iRefSFR = -1
        self.iRefEBmV = -1
        self.iRefSigma = -1
        
    def setMuse2RefFileType(self):
        self.nRefValues = 2
        
        self.iRefID = 0
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
    
    def setPFS201703RefFileType(self):
        """
        for pfs 201703 dataset with avstars and avgas
        """
        self.nRefValues = 10
        self.iRefZ = 1
        self.iRefMag = 2
        self.iRefFlag = -1
        self.iRefSFR = 5
        self.iRefEBmV = 4
        self.iRefSigma = 6
        self.iRefLogHalpha = 7
        self.iRefAvStars = 8
        self.iRefAvGas = 9
        
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
        
    def setSIMUEuclid2017RefFileType(self):
        self.nRefValues = 5
        self.iRefZ = 1
        self.iRefMag = 2
        self.iRefFlag = -1
        self.iRefSFR = 3
        self.iRefEBmV = -1
        self.iRefISM = -1
        self.iRefSigma = -1
        self.iRefLogHalpha = 4    
        
    def setAriadneRefFileType(self):
        self.nRefValues = 6
        self.iRefZ = 1
        self.iRefMag = 2
        self.iRefFlag = -1
        self.iRefSFR = 3
        self.iRefEBmV = -1
        self.iRefISM = -1
        self.iRefSigma = 5
        self.iRefLogHalpha = 4   
        
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
            if not lineStr.startswith('#') and not lineStr=="":
                #print lineStr
                data = lineStr.split("\t")
                if len(data)<2:
                    data = lineStr.split(" ")
                data = [r for r in data if r != '']
                
                _id = -1
                if self.iRefID!=-1:
                    _id = str(data[self.iRefID])
                
                _redshift = -1
                if self.iRefZ!=-1:
                    _redshift = float(data[self.iRefZ])
                
                _mag = -1
                if self.iRefMag!=-1:
                    _mag = float(data[self.iRefMag])
                
                _flag = -1
                if self.iRefFlag!=-1:
                    _flag = str(data[self.iRefFlag])
                
                _sfr = -1
                if self.iRefSFR!=-1:
                    _sfr = float(data[self.iRefSFR])
                
                _ebmv = -1
                if self.iRefEBmV!=-1:
                    _ebmv = float(data[self.iRefEBmV])
                
                _sigma = -1
                if self.iRefSigma!=-1:
                    _sigma = float(data[self.iRefSigma])
                
                _loghalpha = -1
                if self.iRefLogHalpha!=-1:
                    _loghalpha = float(data[self.iRefLogHalpha])
                
                _avstars = -1
                if self.iRefAvStars!=-1:
                    _avstars = float(data[self.iRefAvStars])
                
                _avgas = -1
                if self.iRefAvStars!=-1:
                    _avgas = float(data[self.iRefAvGas])
                
                #filtering
                enableFiltering = False
                if enableFiltering:
                    if 0 and _loghalpha<-17.2:
                        continue
                    
                    if 0 and _sigma<75:
                        continue
                    
                    if 1 and (_redshift<2.0 or _redshift>2.5):
                        continue
                    
                    if 0 and _mag<23.0:
                        continue
                
                #populating 
                self.ids.append(_id)
                self.redshifts.append(_redshift)
                self.mags.append(_mag)
                self.flags.append(_flag)
                self.sfrs.append(_sfr)
                self.ebmvs.append(_ebmv)
                self.sigmas.append(_sigma)
                self.loghalphas.append(_loghalpha)
                self.avstars.append(_avstars)
                self.avgas.append(_avgas)
                
        
        self.n = len(self.ids)
            
    def fakeLoad(self, idList):
        """
        populate ref with fake data
        """ 
        print("INFO: loadref with type={}".format(self.rtype))

        for k in range(len(idList)):
            val = -1
            self.ids.append(idList[k])
            self.redshifts.append(float(val))
            self.mags.append(float(val))
            self.flags.append(str(val))
            self.sfrs.append(float(val))
            self.ebmvs.append(float(val))
            self.sigmas.append(float(val))
            self.loghalphas.append(float(val))
        
        self.n = len(self.ids)
        
        print("INFO: fakeload n={}".format(self.n))
              
    def checkUnicityAllIds(self):
        """
        verify that there are no duplicates when checking inclusion of one id versus another id
        """    
        pb_ids = []
        print("INFO: checking for unicity")          
        for k, _idThis in enumerate(self.ids): 
            if self.findIdx(_idThis)<0:
                print("Inclusion Id problem detected for id={}".format(_idThis))
                if _idThis not in pb_ids:
                    pb_ids.append(_idThis)
        
        if len(pb_ids)<1:
            print("INFO: unicity: no problem found.")
        else:
            for a in pb_ids:
                print("{}\tis_non_unique".format(a))
              
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
        export in the simueuclid2016 format by default, coz it's the most exhaustive format
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
        
    def convertRedshift(self, initialLbda=6562.8, finalLbda=3200.0, idsuffix="c"):
        """
        convert the redshift for lbda center value exchange.
        adds a suffix to the converted ids
        """
        for k, _id in enumerate(self.ids):            
            newRedshift = (1+self.redshifts[k])*initialLbda/finalLbda-1
            self.redshifts[k] = newRedshift
            if not idsuffix=="":
                self.ids[k] = "{}{}".format(self.ids[k], idsuffix)
    
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
            print("ERROR : multiple index found for a given nametag!!".format())
            print("ERROR : idNameTag = {}\n".format(idNameTag))
            print("ERROR : p = {}\n".format(p))
            return -1
        else:
            print("ERROR : index not found for idNameTag={} in the reference ids list!!".format(idNameTag))
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
        
    def plotHistHalphaFlux(self):
        mvect = self.loghalphas
        yvect = [1.0 for a in range(len(mvect))]
        
        print("Plotting versus logfhalpha")
        outputDirectory = os.path.split(self.referencepath)[0]
        outFileNoExt = 'stats_versusLogFHAlpha_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        enablePlot = 0
        nPercentileDepth = 1
        exportType = 'png'
        lstats.PlotAmazedVersusBinsHistogram(yvect, mvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, exportType=exportType, mtype='LOGFHALPHA', nPercentileDepth=nPercentileDepth) 
                
    def plotHistRedshift(self):
        mvect = self.redshifts
        yvect = [1.0 for a in range(len(mvect))]
        
        print("Plotting versus redshift")
        outputDirectory = os.path.split(self.referencepath)[0]
        outFileNoExt = 'stats_versusRedshift_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        enablePlot = 0
        nPercentileDepth = 1
        exportType = 'png'
        lstats.PlotAmazedVersusBinsHistogram(yvect, mvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, exportType=exportType, mtype='REDSHIFT', nPercentileDepth=nPercentileDepth) 
        
    def plotHistMag(self):
        mvect = self.mags
        yvect = [1.0 for a in range(len(mvect))]
        
        print("Plotting versus mag")
        outputDirectory = os.path.split(self.referencepath)[0]
        outFileNoExt = 'stats_versusMag_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        enablePlot = 0
        nPercentileDepth = 1
        exportType = 'png'
        lstats.PlotAmazedVersusBinsHistogram(yvect, mvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, exportType=exportType, mtype='MAG', nPercentileDepth=nPercentileDepth) 
        
    def plotHistSigmaEL(self):
        mvect = self.sigmas
        yvect = [1.0 for a in range(len(mvect))]
        
        print("Plotting versus sigma")
        outputDirectory = os.path.split(self.referencepath)[0]
        outFileNoExt = 'stats_versusSigmaEL_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        enablePlot = 0
        nPercentileDepth = 1
        exportType = 'png'
        lstats.PlotAmazedVersusBinsHistogram(yvect, mvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, exportType=exportType, mtype='SIGMA', nPercentileDepth=nPercentileDepth) 
        
    def plotHistSFR(self):
        mvect = self.sfrs
        yvect = [1.0 for a in range(len(mvect))]
        
        print("Plotting versus sfr")
        outputDirectory = os.path.split(self.referencepath)[0]
        outFileNoExt = 'stats_versusSFR_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        enablePlot = 0
        nPercentileDepth = 1
        exportType = 'png'
        lstats.PlotAmazedVersusBinsHistogram(yvect, mvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, exportType=exportType, mtype='SFR', nPercentileDepth=nPercentileDepth) 
        
    def plotFhaVersusSFR(self):
        xvect = self.sfrs
        yvect = [a for a in self.loghalphas]
        
        plt.figure()
        plt.plot(xvect, yvect, 'x')
        plt.grid()
        plt.xlabel("SFR")
        plt.ylabel("reference {}=f({})".format("logfha", "sfr"))
        plt.xscale('log')
        plt.show() 
        
        
def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--referencePath", dest="referencePath", default="",
                    help="path to the reference file to be loaded")    
                    
    parser.add_argument("-t", "--type", dest="type", default="pfs",
                    help="reference file type, choose between 'simple', 'vvds1', 'vvds2' or 'pfs', 'pfs201703', 'vuds', 'simulm', 'simueuclid2016', 'keck', 'ariadne', 'no'")
          
    options = parser.parse_args()
    #print(options)

    if os.path.exists(options.referencePath) :
        rpath = os.path.abspath(options.referencePath)
        print('using full path: {0}'.format(rpath))
        ref = Reference(rpath, options.type)

        enableConvertAndExport = False
        if enableConvertAndExport:
            print("INFO: Converting redshifts...")
            #finalLbda=3200.0 #foro2
            #id_suffix="_zforo2"
            #finalLbda=5200.0 #for zoff1
            #id_suffix="_zoff1"
            #finalLbda=4150.0 #for zoff2
            #id_suffix="_zoff2"
            finalLbda=6562.8
            #id_suffix="_zforha"
            id_suffix="_zlt0p95"
            
            ref.convertRedshift(initialLbda=6562.8, finalLbda=finalLbda, idsuffix=id_suffix)
            exportFileName = "{}_converted{}.txt".format(os.path.splitext(os.path.basename(rpath))[0], id_suffix)
            exportFilePath = os.path.join(os.path.split(rpath)[0], exportFileName)
            print("INFO: export path is {}".format(exportFilePath))
            if os.path.exists(exportFilePath):
                print("Export path already exists, aborting export...")
            else:
                ref.export(exportFilePath)
        
        ref.checkUnicityAllIds()
        ref.plot()
        #ref.plotHistHalphaFlux()
        #ref.plotHistMag()
        #ref.plotHistRedshift()
        #ref.plotHistSigmaEL()
        #ref.plotHistSFR()
        #ref.plotFhaVersusSFR()
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
    

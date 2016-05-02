# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 16:17:15 2015

@author: aschmitt
"""

import os
from astropy.io import ascii
import numpy as np
import json

import chisquare as chisq

class ResParser(object):
    def __init__(self, respath):
        self.logTagStr = "ResParser"
        self.respath = respath
        self.name = os.path.basename(respath)
        
        ##### results architecture config #####
        # in res directory
        self.redshiftpath = os.path.join(self.respath,"redshift.csv")
        self.configpath = os.path.join(self.respath,"config.txt") 
        self.parameterspath = os.path.join(self.respath,"parameters.json") 
        # template categories
        self.config_tplCategories = ['emission', 'galaxy', 'qso', 'star']
        # in stats directory
        self.statsdirpath = os.path.join(self.respath,"stats")
        self.diffpath = os.path.join(self.statsdirpath,"diff.txt")
        self.failurespath = os.path.join(self.statsdirpath,"failures.txt")

        
        # intermediate directory
        methodParam = self.getParameterVal('method')
        print("using method param : {}".format(methodParam))
        if methodParam=="decisionaltreeb" or methodParam=="amazed0_2":
            self.chisquarerelpath = "dtreeBsolve.chisquare2solve.chisquare.csv" 
            self.chisquarencrelpath = "dtreeBsolve.chisquare2solve.chisquare_nocontinuum.csv"
        else:
            self.chisquarerelpath = "chisquare2solve.chisquare.csv" 
            self.chisquarencrelpath = "chisquare2solve.chisquare_nocontinuum.csv"

        self.chisquarecrelpath = "dtreeBsolve.chisquare2solve.chisquare_continuum.csv"   
        self.corrrelpath = "correlationsolve.correlation.csv"    
        
        
        if methodParam=="linemodelsolve" or methodParam=="linemodel":
            self.linemodelrelpath = "linemodelsolve.linemodel.csv" 
        if methodParam=="linemodeltplshapesolve" or methodParam=="linemodeltplshape":
            self.linemodelrelpath = "linemodeltplshapesolve.linemodel.csv" 
        if methodParam=="decisionaltreeb":
            self.linemodelrelpath = "dtreeBsolve.linemodel.csv" 
                
        self.displaysOutputPath = os.path.join(self.respath,"displays")
           
        
        self.detectedpeakcatalogrelpath = 'linematching2solve.peakdetection.csv'
        self.detectedlinecatalogrelpath = 'linematching2solve.raycatalog.csv'
        self.lineMatchingResultrelpath = 'linematching2solve.raymatching.csv'
        
        self.diffloaded = False
        self.diffData = []
        
        self.failuresloaded = False
        self.failuresData = []

    def __str__(self):
        a = "\nResParser: {0}\n".format(self.name)
        a = a + ("    respath = {0}\n".format(self.respath))
        a = a + ("    redshiftpath = {0}\n".format(self.redshiftpath))
        a = a + ("    configpath = {0}\n".format(self.configpath))
        a = a + ("    diffpath = {0}\n".format(self.diffpath))
        a = a + ("    failurespath = {0}\n".format(self.failurespath))
        a = a + ("\n")
        
        return a  
     
    def getStatsDirPath(self):
        """
        returns the stats path
        """
        return self.statsdirpath
    
    def getDisplaysDirPath(self):
        """
        returns the displays path
        """
        if not os.path.exists(self.displaysOutputPath):
            os.makedirs(self.displaysOutputPath)
        return self.displaysOutputPath
    
    def getDiffPath(self):
        """
        returns the diff path
        """
        return self.diffpath
        
    def loadDiffOld(self):
        """
        load the diff data from file
        """
        fname = self.diffpath
        f = open(fname, 'r')
        dataStr = f.read()
        f.close()
        self.diffData = ascii.read(dataStr, data_start=0, delimiter="\t")
        self.diffloaded = True
        
    def loadDiff(self):
        """
        load the diff data from file
        """ 
        fname = self.diffpath
        f = open(fname)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                data = lineStr.split("\t")
                data = [r for r in data if r != '']
                #print len(data)
                if(len(data) >=9):
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = int(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = str(data[6])
                    d7 = str(data[7])
                    d8 = float(data[8])
                    d = [d0, d1, d2, d3, d4, d5, d6, d7, d8]
                    self.diffData.append(d)
                elif (len(data) ==7):
                    #Spectrum ID	MAGI	ZREF	ZCALC	MERIT	METHOD	DIFF
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = float(data[3])
                    d4 = float(data[4])
                    d5 = str(data[5])
                    d6 = float(data[6])
                    d = [d0, d1, d2, d3, d4, d5, d6]
                    self.diffData.append(d)
                elif (len(data) ==8):
                    #Spectrum ID	MAGI	ZREF	ZFLAG	ZCALC	MERIT	TPL	METHOD	DIFF
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = float(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = str(data[6])
                    d7 = float(data[7])
                    d = [d0, d1, d2, d3, d4, d5, d6, d7]
                    self.diffData.append(d)
        self.diffloaded = True           
        f.close()
        
    def getDiffLine(self, k):
        """
        returns the diff Line by Id
        """
        if not self.diffloaded:
            self.loadDiff()
        return self.diffData[k]
    
    def getDiffSize(self):
        """
        returns the diff data size
        """
        if not self.diffloaded:
            self.loadDiff()
        return len(self.diffData)
    
    def getConfigVal(self, tag):
        """
        returns the string value of the field corresponding to the tag in the config file 
        """
        strVal = ""
        filename = self.configpath
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            data = lineStr.split("=")
            if(len(data) >=2):
                if(tag == data[0]):
                    strVal = data[1]
                    break
        f.close()
        return strVal
    
    def getParameterVal(self, tag1, tag2="", tag3="", tag4="", tag5=""):
        """
        returns the string value of the field corresponding to the tag in the parameter file 
        """
        strVal = ""
        filename = self.parameterspath
        with open(filename) as data_file:    
            data = json.load(data_file)
        
        if tag2=="":
            strVal = data[tag1]
        elif tag3=="":
            strVal = data[tag1][tag2]
        elif tag4=="":
            strVal = data[tag1][tag2][tag3]
        elif tag5=="":
            strVal = data[tag1][tag2][tag3][tag4]
        else:
            strVal = data[tag1][tag2][tag3][tag4][tag5]
        return strVal
        
    def getSpectrumlistline(self, tag):
        """
        returns the string array (1xi, with i in {1,2}) of the field corresponding to the tag in the spectrumlist file 
        """
        outVal = ""
        filename = self.getConfigVal('input')
        f = open(filename)
        #print filename
        for line in f:
            if line.startswith("#"):
                continue
            if line == "":
                continue
            lineStr = line.strip()
            #print lineStr
            data = lineStr.split("\t")
            data = [r for r in data if r != '']
            if(len(data) >=1):
                if not (data[0].find(tag)==-1):
                    outVal = data
                    break
        f.close()
        return outVal
  
    def getRedshiftVal(self, spcnametag):
        """
        returns the redshift value for the spectrum corresponding to the nametag in the output redshift file 
        """
        spcnametag = os.path.splitext(spcnametag)[0]
        #print("getRedshiftVal : spcnametag = {}".format(spcnametag))
        floatVal = ""
        filename = self.redshiftpath
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            data = lineStr.split("\t")
            if(len(data) >=4):
                if(spcnametag == data[0]):
                    floatVal = float(data[1])
                    break
        f.close()
        return floatVal
        
    def getRedshiftTpl(self, spcnametag):
        """
        returns the template selected for the spectrum corresponding to the nametag in the output redshift file 
        """
        #override:
        if 0:
            tplnametag = "zcosmos_red.txt"
            tplnametag = "NEW_Sbc_extended.dat"
            tplnametag = "COMBINE-ave-Lya-abs-AND-StarBurst1.txt"
            tplnametag = "Rebinned-NEW-E-extendeddataExtensionData"
            return tplnametag     
            
        strVal = ""
        filename = self.redshiftpath
        #print("getRedshiftTpl, filename : {}".format(filename))
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            data = lineStr.split("\t")
            if(len(data) >=4):
                #print("getRedshiftTpl, data[0] : {}".format(data[0]))
                if(spcnametag in data[0]):
                    strVal = data[3]
                    break
        f.close()
        return strVal
        
    def getAutoTplFullPath(self, spcnametag, idxExtrema=0):
        """
        """ 
        if os.path.splitext(spcnametag)[1].lower()==".fits":
            spcnametag = os.path.splitext(spcnametag)[0]
        method = self.getConfigVal('method')
        print("method found in config is: {}".format(method))
        path = os.path.join(self.respath, spcnametag)

        tplpath = ""
        if method == "linemodel":
            name = "linemodelsolve.linemodel_spc_extrema_{}.csv".format(idxExtrema)
            tplpath = os.path.join(path,name)
        if method == "linemodeltplshape":
            tplnametag = self.getRedshiftTpl(spcnametag)
            pathTplChi = os.path.join(path, tplnametag)
            name = "linemodeltplshapesolve.linemodel_spc_extrema_{}.csv".format(idxExtrema) 
            tplpath = os.path.join(pathTplChi,name)
            #tplpath = ""
        elif method == "decisionaltreeb" or method.lower() == "amazed0_2":
            name = "dtreeBsolve.linemodel_spc_extrema_{}.csv".format(idxExtrema)
            tplpath = os.path.join(path,name)
            #tplpath = ""
            
        return tplpath
        
    def getTplFullPath(self, tplnametag):
        """
        """
        tplnametag_noext = self.getWithoutExt(tplnametag)
        strVal = ""
        tplRootPath = self.getConfigVal('templatedir')
        for cat in self.config_tplCategories:
            tplPath = os.path.join(tplRootPath, cat)
            if os.path.exists(tplPath):
                for file in sorted(os.listdir(tplPath)):
                    if tplnametag_noext == self.getWithoutExt(file):
                        strVal = os.path.join(tplPath, file)
                        break
                    
        return strVal
        
    def getSpcFullPath(self, spcnametag):
        """
        """
        #print("spcnametag: {}".format(spcnametag))
        spcnametag_noext = self.getWithoutExt(spcnametag)
        #print("spcnametag_noext: {}".format(spcnametag_noext))
        strVal = ""
        spcPath = self.getConfigVal('spectrumdir')
        if os.path.exists(spcPath):
            for file in sorted(os.listdir(spcPath)):
                #print("file: {}".format(file))
                #print("self.getWithoutExt(file): {}".format(self.getWithoutExt(file)))
                if spcnametag_noext == self.getWithoutExt(file):
                    strVal = os.path.join(spcPath, file)
                    break
                    
        return strVal
        
    def getNoiseFullPath(self, spcnametag):
        """
        """
        #print("getNoiseFullPath: spcnametag: {}".format(spcnametag))
        spectrumlistline = self.getSpectrumlistline(spcnametag)
        print("getNoiseFullPath: spcnametag: {}".format(spectrumlistline))
        if len(spectrumlistline) == 0:
            print("WARNING: unable to find spectrumlistline for the tag: {}".format(spcnametag))
            strVal = ""
        if len(spectrumlistline) == 1:
            strVal = ""
        if len(spectrumlistline) == 2:
            noisetag = spectrumlistline[1]
            noisetag_noext = self.getWithoutExt(noisetag)
            spcPath = self.getConfigVal('spectrumdir')
            if os.path.exists(spcPath):
                for file in sorted(os.listdir(spcPath)):
                    #print("file: {}".format(file))
                    #print("self.getWithoutExt(file): {}".format(self.getWithoutExt(file)))
                    if noisetag_noext == self.getWithoutExt(file):
                        strVal = os.path.join(spcPath, file)
                        break
        
        return strVal
                
    def getCatalogFullPath(self):
        """
        """
        strVal = self.getConfigVal('linecatalog')
        return strVal
    
    
    
    def getWithoutExt(self, tag):
        """
        """
        extensions = ['.fits', '.FITS', '.txt', '.TXT', '.dat', '.DAT']
        for e in extensions:
            if tag.endswith(e):
                iend = len(e)
                data = tag[:-iend]
                return data
        return tag
            
    def getAutoChi2FullPath(self, spcnametag):
        """
        """ 
        if os.path.splitext(spcnametag)[1].lower()==".fits":
            spcnametag = os.path.splitext(spcnametag)[0]
        method = self.getConfigVal('method')
        print("method found in config is: {}".format(method))
        path = os.path.join(self.respath, spcnametag)
        #, chi2type="raw"
        #tplnametag
        chipath = []
        chiname = []
            
        if method == "chisquaresolve" or method == "chisquare2solve":
            tplnametag = self.getRedshiftTpl(spcnametag)
            pathTplChi = os.path.join(path, tplnametag)
            spcComponent = self.getParameterVal('chisquare2solve', 'spectrum', 'component')
            print('component parameter found = {}'.format(spcComponent))
            if spcComponent=="nocontinuum":
                name = "chisquare2solve.chisquare_nocontinuum.csv"
            if spcComponent=="continuum":
                name = "chisquare2solve.chisquare_continuum.csv"
            if spcComponent=="raw":
                name = "chisquare2solve.chisquare.csv"
            #name = "chisquare2solve.chisquare.csv"
            #name = "correlationsolve.correlation.csv"
            #name = "chisquare2solve.chisquare_nocontinuum.csv"
            chipath.append(os.path.join(pathTplChi,name))
            chiname.append(name)
            
        elif method == "correlationsolve":
            tplnametag = self.getRedshiftTpl(spcnametag)
            pathTplChi = os.path.join(path, tplnametag)
            name = "correlationsolve.correlation.csv"
            chipath.append(os.path.join(pathTplChi,name))
            chiname.append(name)
            
        elif method == "decisionaltree7":
            tplnametag = self.getRedshiftTpl(spcnametag)
            pathTplChi = os.path.join(path, tplnametag)
            name = "dtree7solve.correlationsolve.correlation.csv"
            #name = "correlationsolve.correlation.csv"
            #name = "chisquare2solve.chisquare_nocontinuum.csv"
            chipath.append(os.path.join(pathTplChi,name))
            chiname.append(name)
            
        elif method == "blindsolve":
            tplnametag = self.getRedshiftTpl(spcnametag)
            pathTplChi = os.path.join(path, tplnametag)
            name = "blindsolve.correlation.csv"
            #name = "correlationsolve.correlation.csv"
            #name = "chisquare2solve.chisquare_nocontinuum.csv"
            chipath.append(os.path.join(pathTplChi,name))
            chiname.append(name)
            
        elif method == "linemodel":
            name = "linemodelsolve.linemodel.csv"
            #name = "dtreeBsolve.linemodel.csv"    
            chipath.append(os.path.join(path,name))
            chiname.append(name)
            
        elif method == "linemodeltplshape":
            tplnametag = self.getRedshiftTpl(spcnametag)
            pathTplChi = os.path.join(path, tplnametag)
            
            name = "linemodeltplshapesolve.linemodel.csv"  
            chipath.append(os.path.join(pathTplChi,name))
            chiname.append(name)
            
        elif method == "decisionaltreeb" or method.lower() == "amazed0_2":
            
            if 1:
                #name = "dtreebsolve.linemodel.csv"
                name = "dtreeBsolve.resultdtreeBCombined.csv"    
                chipath.append(os.path.join(path,name))
                chiname.append(name)
                
           
            if 1:
                tplnametag = "EW_SB2extended.dat"
                tplnametag = "NEW_Sbc_extended.dat"
                #tplnametag = "NEW_Im_extended_blue.dat"NEW_Sbc_extended
                tplnametag = self.getChi2BestTpl(spcnametag, chi2type="continuum")
                
                print("tpl name tag found = {}".format(tplnametag))
                pathTplChi = os.path.join(path, tplnametag)
                name = "dtreeBsolve.chisquare2solve.chisquare_continuum.csv"
                chipath.append(os.path.join(pathTplChi,name))
                chiname.append(name)            
            if 1:
                tplnametag = self.getChi2BestTpl(spcnametag, chi2type="nocontinuum")
                
                print("tpl name tag found = {}".format(tplnametag))
                pathTplChi = os.path.join(path, tplnametag)
                name = "dtreeBsolve.chisquare2solve.chisquare_nocontinuum.csv"
                chipath.append(os.path.join(pathTplChi,name))
                chiname.append(name)
                
            if 1:
                #name = "dtreebsolve.linemodel.csv"
                name = "dtreeBsolve.linemodel.csv"    
                chipath.append(os.path.join(path,name))
                chiname.append(name)

        return chipath, chiname

         
    def getChi2FullPath(self, spcnametag, tplnametag, chi2type="raw"):
        """
        """   
        dirpath = os.path.join(self.respath, spcnametag)
        #print("dirpath = {0}".format(dirpath))
        if not chi2type=='linemodel':
            dirpath = os.path.join(dirpath, tplnametag)
        #print("dirpath = {0}".format(dirpath))
            
        frelpath = self.chisquarerelpath
        if chi2type=='raw':
            frelpath = self.chisquarerelpath
        elif chi2type=='continuum':
            frelpath = self.chisquarecrelpath
        elif chi2type=='nocontinuum':
            frelpath = self.chisquarencrelpath
        elif chi2type=='corr':
            frelpath = self.corrrelpath 
        elif chi2type=='linemodel':
            frelpath = self.linemodelrelpath
        elif chi2type=='linemodeltplshape':
            frelpath = self.linemodelrelpath 
        
        chi2path = os.path.join(dirpath, frelpath)  
        return chi2path

        
    def getChi2Val(self, spcnametag, zcalc, tplName="", chi2type="raw"):
        """
        returns the chi2 value for the spectrum corresponding to the nametag 
        in the output intermediate chi2 file 
        """
        floatVal = -1  
        if tplName=="":
            tplName = self.getRedshiftTpl(spcnametag);
        filename = self.getChi2FullPath(spcnametag, tplName, chi2type)
        #print('getChi2Val using full path: {0}'.format(filename))
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                data = lineStr.split("\t")
                if(len(data) >=2):
                    if(zcalc -float( data[0]) < 1e-8):
                        floatVal = float(data[1])
                        break
        f.close()
        return floatVal
 
    def getChi2MinValList(self, spcnametag, chi2type="raw", dontloadThres=-1):
        """
        returns the chi2 min value (minimimum over the tpl) list for the spectrum corresponding to the nametag 
        in the output intermediate chi2 file 
        """
        chi2list = []
        if chi2type=='linemodel':
            tplList = ["unusedtplpathlist"]
        else:
            tplList = self.getTplFullPathList()
        ntpl = len(tplList)
        #get the number of redshifts values: nz
        
        nz=0
        xvect=[]
        for x in range(ntpl):
            tplName= os.path.basename(tplList[x])
            chipath = self.getChi2FullPath(spcnametag, tplName, chi2type)
            chi = chisq.ResultChisquare(chipath, stype='undef', dontloadThres=dontloadThres)
            
            print("getChi2MinValList - chi, xvect n = {}".format(len(chi.xvect)))
            if(nz<len(chi.xvect)):
                xvect = chi.xvect
                nz = len(xvect)
        print("getChi2MinValList - nz = {}".format(nz))
        #initialize the array
        fullchi2s = np.ones((nz,ntpl))
        
        for x in range(ntpl):
            tplName= os.path.basename(tplList[x])
            chipath = self.getChi2FullPath(spcnametag, tplName, chi2type)
            chi = chisq.ResultChisquare(chipath, stype='undef', dontloadThres=dontloadThres)
            ymin = 0
            for y in range(nz):
                if xvect[ymin] < chi.xvect[0]:
                    fullchi2s[y,x] = 1e12
                    ymin = ymin+1
                if xvect[y] > chi.xvect[len(chi.xvect)-1]:
                    fullchi2s[y,x] = 1e12
                else:
                    fullchi2s[y,x] = chi.yvect[y-ymin]
        for y in range(nz):
            chi2list.append(min(fullchi2s[y, :]))
            
        return xvect, chi2list       
 
    def getChi2BestTpl(self, spcnametag, chi2type="raw"):
        """
        returns the chi2 best template for the spectrum corresponding to the nametag 
        in the output intermediate chi2 file 
        """
        if chi2type=='linemodel':
            tplList = ["unusedtplpathlist"]
        else:
            tplList = self.getTplFullPathList()
        ntpl = len(tplList)
        #get the number of redshifts values: nz
        
        chi2MinVal=1e12
        bestTplName = ""
        for x in range(ntpl):
            tplName= os.path.basename(tplList[x])
            chipath = self.getChi2FullPath(spcnametag, tplName, chi2type)
            chi = chisq.ResultChisquare(chipath)
            for y in range(len(chi.yvect)):
                if chi2MinVal > chi.yvect[y]:
                    chi2MinVal = chi.yvect[y]
                    bestTplName = tplName
            
        return bestTplName
        
    def getChi2List(self, spcnametag, zcalc, chi2type="raw" ):
        """
        returns the chi2 value for all the templates (list), for the spectrum 
        corresponding to the nametag in the output intermediate chi2 file 
        """
        chi2list = []
        tplList = self.getTplFullPathList()
        ntpl = len(tplList)
        for x in range(ntpl):
            tplName= os.path.basename(tplList[x])
            chi2 = self.getChi2Val(spcnametag, zcalc, tplName, chi2type)
            chi2list.append(chi2)
        return chi2list
        
    def getTplFullPathList(self):
        """
        returns the list of templates full paths
        """
        tplPathList = []
        tplRootPath = self.getConfigVal('templatedir')
        for cat in self.config_tplCategories:
            tplPath = os.path.join(tplRootPath, cat)
            if os.path.exists(tplPath):
                for file in sorted(os.listdir(tplPath)):
                    tplPathList.append(os.path.join(tplPath, file)) 
        return tplPathList
        
    def getDetectedLineCatalogPath(self, spcnametag):
        """
        """   
        dirpath = os.path.join(self.respath, spcnametag)
        #print("dirpath = {0}".format(dirpath))
        dctlgfullpath = os.path.join(dirpath, self.detectedlinecatalogrelpath)
 
        return dctlgfullpath    
        
    def getDetectedPeakCatalogPath(self, spcnametag):
        """
        """   
        dirpath = os.path.join(self.respath, spcnametag)
        #print("dirpath = {0}".format(dirpath))
        peakctlgfullpath = os.path.join(dirpath, self.detectedpeakcatalogrelpath)
 
        return peakctlgfullpath
        
    def getLineMatchingResultPath(self, spcnametag):
        """
        """   
        dirpath = os.path.join(self.respath, spcnametag)
        #print("dirpath = {0}".format(dirpath))
        dctlgfullpath = os.path.join(dirpath, self.lineMatchingResultrelpath)
 
        return dctlgfullpath


if __name__ == '__main__':
    
    rpath = "/home/aschmitt/data/pfs/pfs_lbg/amazed/res_20150706_chisquare"
    print('using full path: {0}'.format(rpath))
    s = ResParser(rpath)
    print(s) 
    print(s.getSpcFullPath('EZ_fits-W-F_3'))
    print(s.getRedshiftVal('EZ_fits-W-F_3'))
    print(s.getRedshiftTpl('EZ_fits-W-F_3'))
    print(s.getTplFullPath(s.getRedshiftTpl('EZ_fits-W-F_3')))
    
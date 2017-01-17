# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 10:49:18 2015

@author: aschmitt
"""
import sys
import os
import inspect
import optparse
import numpy as np
import math
import time
import traceback

import matplotlib as mpl
mpl.use('Qt5Agg')

#mpl.use('Agg') #disable showing the mpl windows
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D

from scipy.interpolate import griddata
from scipy import ndimage
import scipy.optimize as opt
from scipy import interpolate

subfolder = "../stats"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],subfolder)))
if cmd_subfolder not in sys.path:
    #print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
    sys.path.insert(0, cmd_subfolder)
     
import lstats

import resparser as rp
import chisquare as chisq  
import detectedpeakcatalog as dpeakctlg
import detectedcatalog as dctlg 
import matchingsolution as matchsol 
import spectrum as sp
import modelresult

class Result(object):
    def __init__(self, name, zref, zcalc, zdiff, chi2, chi2nc, chi2PerTplZcalc, chi2PerTplZref, chi2ncPerTplZref, refValues):
        self.logTagStr = "Result"
        self.name = name
        self.zcalc = zcalc;
        self.zref = zref;
        self.zdiff = zdiff;
        self.chi2 = chi2;
        self.chi2nc = chi2nc;
        
        self.chi2PerTplZcalc = chi2PerTplZcalc; #list of chi2 val, per template, at z = zref
        self.chi2PerTplZref = chi2PerTplZref; #list of chi2 val, per template, at z = zref
        self.chi2ncPerTplZref = chi2ncPerTplZref; #list of chi2nocontinuum val, per template, at z = zref
        if len(chi2PerTplZref)>0:
            self.tplMissingRate = self.ComputeTplMissingRate(); # rate: percentage of templates for which the chi2 val is not available at z ref, could be due to insufficient overlap

        self.refValues = refValues        
        
    def ComputeTplMissingRate(self):
        ntpl = len(self.chi2PerTplZref)
        if ntpl <1:
            return -1
        nMissing = 0
        for a in range(ntpl):
            if self.chi2PerTplZref[a] > 1e20:
                nMissing +=1
        rate = nMissing/float(ntpl);
        return rate;
                
class ResultList(object):
    def __init__(self, dir, diffthreshold=-1, opt='full', spcName="", methodName="", zrefmin=-1, zrefmax=20, magrefmin=-100.0, magrefmax=100, sfrrefmin=-1.0, sfrrefmax=1e4):
        self.logTagStr = "ResultList"
        self.dir = dir
        self.name = os.path.basename(self.dir)
        self.diffthreshold = diffthreshold  
        self.zrefmin = zrefmin
        self.zrefmax = zrefmax
        self.magrefmin = magrefmin
        self.magrefmax = magrefmax
        self.sfrrefmin = sfrrefmin
        self.sfrrefmax = sfrrefmax
        
        self.opt = opt
        self.statsdir = ""
        self.list = [] 
        self.n = -1;

        self.analysisoutputdir = ""
        
        
        self.load(spcName=spcName, methodName=methodName)
    
    def load(self, spcName="", methodName=""):
        self.resParser = rp.ResParser(self.dir)
        print(self.resParser) 
        self.statsdir = self.resParser.getStatsDirPath()
        self.analysisoutputdir = os.path.join(self.statsdir, "analysis")
        diffthresStr = "diffthres_{}".format(self.diffthreshold)
        if(self.diffthreshold<0):
            diffthresStr = "diffthres_none"
        self.analysisoutputdir = os.path.join(self.analysisoutputdir, diffthresStr)
        if not os.path.exists(self.analysisoutputdir):
            os.makedirs(self.analysisoutputdir)
    
        n = self.resParser.getDiffSize() 
        #n = 4
        print('Number of spectra is: {0}'.format(n))
    
        self.n = 0;
        for x in range(0, n):
            line = self.resParser.getDiffLine(x)
            name = line[0] 
            zref = line[2]

            magref = line[1]
            sfrref = line[9] 
            ebmvref = line[10]            
            
            zcalc = line[4]
            zdiff = line[4]-line[2]
            #print('zcalc is: {0}'.format(zcalc)) 
            method = line[6]
            
            elvelocity = line[11]
            
            enableShowDetails = False
            if enableShowDetails:
                print('name is: {0}'.format(name))
                print('zref is: {0}'.format(zref))
                print('zcalc is: {0}'.format(zcalc)) 
                print('self.diffthreshold is: {0}'.format(self.diffthreshold))
                print('method is: {0}'.format(method)) 
                print('elvelocity is: {0}'.format(elvelocity)) 
            if self.opt=='full':
                chi2 = self.resParser.getChi2Val(name, zcalc, "", chi2type="raw")
                chi2nc = self.resParser.getChi2Val(name, zcalc, "", chi2type="nocontinuum")
                chi2PerTplZcalc = self.resParser.getChi2List(name, zcalc, chi2type="raw")
                chi2PerTplZref = self.resParser.getChi2List(name, zref, chi2type="raw")
                chi2ncPerTplZref = [];#self.resParser.getChi2List(name, zref, chi2type="nocontinuum")
            else:
                chi2 = -1
                chi2nc = -1
                chi2PerTplZcalc = []
                chi2PerTplZref = []
                chi2ncPerTplZref = []
            
            accepted = True
            if not (spcName=="" or spcName==name): 
                accepted = False
                if enableShowDetails:
                    print("Rejected by name")
            if not (methodName=="" or methodName==method):  
                accepted = False
                if enableShowDetails:
                    print("Rejected by methodName")
            if not abs(zdiff)>=self.diffthreshold:
                accepted = False
                if enableShowDetails:
                    print("Rejected by diffthreshold")
            if not (zref>=self.zrefmin and zref<=self.zrefmax):
                accepted = False
                if enableShowDetails:
                    print("Rejected by zref")
            if not (magref>=self.magrefmin and magref<=self.magrefmax):
                accepted = False
                if enableShowDetails:
                    print("Rejected by magref")
            if not (sfrref>=self.sfrrefmin and sfrref<=self.sfrrefmax):
                accepted = False
                if enableShowDetails:
                    print("Rejected by sfrref")
                
            if accepted:
                refValues = {'elvelocity': elvelocity, 'mag': magref, 'sfr': sfrref, 'ebmv': ebmvref}
                res = Result(name, zref, zcalc, zdiff, chi2, chi2nc, chi2PerTplZcalc, chi2PerTplZref, chi2ncPerTplZref, refValues)
                self.list.append(res)
                self.n +=1
            #print('{0}/{1} results loaded: last chi2 = {2}'.format(self.n,x, chi2))
            if enableShowDetails:
                print('{0}/{1} results loaded\n'.format(self.n,x+1))
                
        print('{0} results loaded.'.format(self.n))
        
    def getZCandidatesFromDiff(self, indice=0):
        return [self.list[indice].zref, self.list[indice].zcalc]
    
    def getZCandidatesFromChi2Extrema(self, indice=0, chi2Type = "raw", nextrema=5):
        redshiftslist = []
        spcName = self.list[indice].name
        s = rp.ResParser(self.dir)
        tplpaths = s.getTplFullPathList()
        for a in range(len(tplpaths)):
            filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), chi2Type)
            chi2 = chisq.ResultChisquare(filepath)
            redshiftstpl = chi2.getFluxExtrema2(nextrema,  enablePlot=False)
            for b in range(nextrema):
                redshiftslist.append(redshiftstpl[b])
        redshiftsUnsorted = list(set(redshiftslist))
        sortId=np.argsort(redshiftsUnsorted)
        redshifts = [redshiftsUnsorted[b] for b in sortId]
        return redshifts
    
    def getZCandidatesFromAmazedChi2Extrema(self, indice=0, chi2Type = "raw", nextrema = 5, enableZrangePerTpl=False):
        """
        retrieve the amazed intermediate results redshift candidates
        - nextrema is the number of candidates to retrieve by order of merit
        - for template fitting [chi2, linemodel-tplshape] (there are n merit functions to consider), the nextrema first 
            extrema are extracted by merit for each templates and ordered in a 
            global list, which is then cropped to contain only nextrema values
            todo: maybe the extrema duplicates should be removed, eg. if abs(zTpl1-zTpl2)<threshold tbd
        """
        redshiftslist = []
        meritslist = []        
        spcName = self.list[indice].name
        print("spcname = {}".format(spcName))
        s = rp.ResParser(self.dir)
        if not chi2Type=='linemodel':
            tplpaths = s.getTplFullPathList()
        else:
            tplpaths = [''];
        if len(tplpaths)<1:
            print("Problem while getting the templates directory path...")

        nextremaPerTpl = max(1,nextrema)
        for a in range(len(tplpaths)):
            #filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), "nocontinuum")
            #filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), "raw")
            filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), chi2Type)
            if not os.path.exists(filepath):
                print("Problem while retrieving chi2 filepath.. using: {}".format(filepath))
                continue
            else:
                print("using Chi2 file path : ".format(filepath))
            chi2 = chisq.ResultChisquare(filepath)
            
            if not enableZrangePerTpl:
                redshiftstpl = chi2.amazed_extrema
                
                if len(redshiftstpl)<nextremaPerTpl:
                    print("number of amazed_extrema found is too small = n={}, vs nextrema={}".format(len(redshiftstpl), nextremaPerTpl))
                    return [-1],[-1]
            else:
                tpltag = os.path.basename(tplpaths[a])
                zrange = self.getZrangeForTemplate(tpltag)
                redshiftstpl = [a for a in chi2.amazed_extrema if a>=zrange[0] and a<= zrange[1]]
                if len(redshiftstpl)<nextremaPerTpl:
                    print("number of amazed_extrema found is too small = n={}, vs nextrema={}".format(len(redshiftstpl), nextremaPerTpl))
                    continue
            print("amazed_extrema found = {}".format(redshiftstpl))
            
            
            meritstpl = []
            npx = np.copy(chi2.xvect)
            for b in range(nextremaPerTpl):
                xfind = np.abs(npx-redshiftstpl[b])
                ind = np.argmin(xfind)
                meritstpl.append(chi2.yvect[ind])
                #print("merit for z={} : {}".format(redshiftstpl[b], chi2.yvect[ind]))
            for b in range(nextremaPerTpl):
                redshiftslist.append(redshiftstpl[b])
                meritslist.append(meritstpl[b])
        if nextrema > 0:
            redshiftsUnsorted = list((redshiftslist))
            meritsUnsorted = list((meritslist))
            sortId=np.argsort(meritsUnsorted)
            if chi2Type == "corr":
                sortId = sortId[::-1]
            redshifts_all = [redshiftsUnsorted[b] for b in sortId]
            merits_all = [meritsUnsorted[b] for b in sortId]

            if False and len(redshifts_all)>nextrema: #disabled the final list crop to nextrema 
                redshifts = [a for i, a in enumerate(redshifts_all) if i < nextrema]
                merits = [a for i, a in enumerate(merits_all) if i < nextrema]
            else:
                redshifts = redshifts_all
                merits = merits_all
            #print("getZCandidatesFromAmazedChi2Extrema, len(redshifts={}".format(len(redshifts))) 
            
            
        elif nextrema == 0:
            redshiftsUnsorted = list((redshiftslist))
            meritsUnsorted = list((meritslist))
            print("meritsUnsorted = {}".format(meritsUnsorted))
            sortId=np.argsort(meritsUnsorted)
            if chi2Type == "corr":
                redshifts = [redshiftsUnsorted[sortId[len(sortId)-1]]]
                merits = [meritsUnsorted[sortId[len(sortId)-1]]]
            else:
                redshifts = [redshiftsUnsorted[sortId[0]]]
                merits = [meritsUnsorted[sortId[0]]]
        
        return redshifts, merits

    def getBestZinZrange(self, indice=0, chi2Type = "raw", enableZrangeFilter=False): 
        bestz = -1
        bestmerit = 1e12
        besttpl = "-1"
        
        spcName = self.list[indice].name
        print("spcname = {}".format(spcName))
        s = rp.ResParser(self.dir)
        if not chi2Type=='linemodel':
            tplpaths = s.getTplFullPathList()
        else:
            tplpaths = [''];
        if len(tplpaths)<1:
            print("Problem while getting the templates directory path...")

        for a in range(len(tplpaths)):
            filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), chi2Type)
            if not os.path.exists(filepath):
                print("Problem while retrieving chi2 filepath.. using: {}".format(filepath))
                continue
            else:
                print("using Chi2 file path : ".format(filepath))
            chi2 = chisq.ResultChisquare(filepath)
            
            tpltag = os.path.basename(tplpaths[a])
            zrange = [-1000.0, 1000.0] #using full zrange by default
            if enableZrangeFilter:
                zrange = self.getZrangeForTemplate(tpltag)
            print("using zrange = {}".format(zrange))
            z, merit = chi2.getBestZ_byFluxMin(zrange[0], zrange[1])
            print("found best redshift for this tpl = {}, with merit = {}".format(z, merit))
            
            if bestmerit>merit:
                bestz=z
                bestmerit=merit
                besttpl = tpltag
                
        print("best z found = {}, with merit = {}".format(bestz, bestmerit))
        return bestz, bestmerit, besttpl
        
    def to_flag(self):
        thres = 0.01
        flags = []
        for x in range(0,self.n):
            if abs(self.list[x].zdiff)>thres:
                flags.append(0)
            else:
                flags.append(1)
        return flags
        
    def to_color(self):
        thres = 0.01
        colors = []
        for x in range(0,self.n):
            if abs(self.list[x].zdiff)>thres:
                colors.append('r')
            else:
                colors.append('k')
        return colors
        
    def chi2_to_liste(self):
        chi2 = [a.chi2 for a in self.list]
        chi2nc = [a.chi2nc for a in self.list]
        return [chi2, chi2nc]

    def chi2_export_liste(self):
        fname = os.path.join(self.statsdir, "resultlist.csv")
        f = open(fname, 'w')
        for item in self.to_liste():
            f.write("%s\n" % item)
        f.close()
        
    def chi2_plot(self):
        Y = self.chi2_to_liste()
        c = self.to_color()
        #print len(Y)
        #print len(Y[0])
        xplot = Y[0];# - np.mean(Y[0],axis=0)
        yplot = Y[1];# - np.mean(Y[1],axis=0)
        for x in range(0,self.n):
            plt.plot(xplot[x],yplot[x], 'x', color=c[x])
        plt.show()
    
    def tplMissingRate_to_liste(self):
        tplMissingRateLst = [a.tplMissingRate for a in self.list]
        return tplMissingRateLst
        
    def tplMissingRate_plot(self):
        yvect = self.tplMissingRate_to_liste()
        c = self.to_color()
        
        for x in range(0,self.n):
            plt.plot(x,yvect[x], 'x', color=c[x])
        plt.show()
        
        fig = plt.figure('Template Missing Rate Histogram')
        ax = fig.add_subplot(111)
        #ax.plot(vectErrorBins, yVectErrorBins)
        plt.hist(yvect, len(self.list[0].chi2PerTplZref), normed=0, histtype='stepfilled', cumulative=False)
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        #ax.set_xscale('log')
        plt.grid(True) # Affiche la grille
        #pp.legend(('cos','sin'), 'upper right', shadow = True)
        plt.ylabel('Count')
        plt.xlabel('Missing rate')
        plt.show()

    def getBestTplZRef(self):
        bestTpl = []
        s = rp.ResParser(self.dir)
        #print(s) 
        tplPaths = s.getTplFullPathList()
        tplNames = [os.path.basename(a) for a in tplPaths]  
        for x in range(0,self.n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            imin = np.argmin(self.list[x].chi2PerTplZref)
            bestTpl.append(tplNames[imin])
        return bestTpl
    
    def getBestChi2ZRef(self):
        bestChi2 = []
        s = rp.ResParser(self.dir)
        #print(s) 
        tplPaths = s.getTplFullPathList()
        tplNames = [os.path.basename(a) for a in tplPaths]  
        for x in range(0,self.n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            imin = np.argmin(self.list[x].chi2PerTplZref)
            bestChi2.append(self.list[x].chi2PerTplZref[imin])
        return bestChi2
        
    def getBestTplZCalc(self):
        bestTpl = []
        s = rp.ResParser(self.dir)
        #print(s) 
        tplPaths = s.getTplFullPathList()
        tplNames = [os.path.basename(a) for a in tplPaths]  
        for x in range(0,self.n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            imin = np.argmin(self.list[x].chi2PerTplZcalc)
            bestTpl.append(tplNames[imin])
        return bestTpl 

        
    def getZrangeForTemplate(self, tplTag):
        zrangePerTemplate = {}
        
        if 0:
            #DEPRECATED: ExtendedGalaxyEL2/, for lambda range = [4500, 9400]
            zrangePerTemplate["BulgedataExtensionData.dat"]=[0.0, 2.076]
            zrangePerTemplate["EdataExtensionData.dat"]=[0.0, 0.581]
            zrangePerTemplate["EllipticaldataExtensionData.dat"]=[0.0, 2.076]
            zrangePerTemplate["E_reddataExtensionData.dat"]=[0.0, 2.081]
            zrangePerTemplate["EW_SB2extended.dat"]=[0.261, 3.607]
            zrangePerTemplate["NEW_E_extendeddataExtensionData.dat"]=[0.0, 4.516]
            zrangePerTemplate["NEW_Im_extended.dat"]=[0.261, 3.146]
            zrangePerTemplate["NEW_Im_extended_blue.dat"]=[0.261, 2.568]
            zrangePerTemplate["NEW_Sbc_extended.dat"]=[0.261, 3.146]
            zrangePerTemplate["s0dataExtensionData.dat"]=[0.0, 2.076]
            zrangePerTemplate["sadataExtensionData.dat"]=[0.0, 2.076]
            zrangePerTemplate["sbdataExtensionData.dat"]=[0.0, 2.076]
            zrangePerTemplate["Scd.txt"]=[0.261, 0.583]
            zrangePerTemplate["StarBurst1.txt"]=[0.261, 0.583]
            zrangePerTemplate["StarBurst2.txt"]=[0.261, 0.583]
            zrangePerTemplate["StarBurst3.txt"]=[0.261, 0.583]
            zrangePerTemplate["vvds_reddestdataExtensionData.dat"]=[0.0, 0.696]
            zrangePerTemplate["zcosmos_red.txt"]=[2.34, 3.48]
        
        if 0:
            #from Amazed outputs, chi2nc 20160321: ExtendedGalaxyEL2/, for lambda range = [3600, 9460]
            zrangePerTemplate["BulgedataExtensionData.dat"]=[0.0, 1.9163]
            zrangePerTemplate["EdataExtensionData.dat"]=[0.0, 0.4989]
            zrangePerTemplate["EllipticaldataExtensionData.dat"]=[0.0, 1.9164]
            zrangePerTemplate["E_reddataExtensionData.dat"]=[0.0, 1.9212]
            zrangePerTemplate["EW_SB2extended.dat"]=[0.0, 3.3668]
            zrangePerTemplate["NEW_E_extendeddataExtensionData.dat"]=[0.0, 4.2290]
            zrangePerTemplate["NEW_Im_extended.dat"]=[0.0, 2.9304]
            zrangePerTemplate["NEW_Im_extended_blue.dat"]=[0.0, 2.3820]
            zrangePerTemplate["NEW_Sbc_extended.dat"]=[0.0, 3.2469]
            zrangePerTemplate["s0dataExtensionData.dat"]=[0.0, 1.9164]
            zrangePerTemplate["sadataExtensionData.dat"]=[0.0, 1.9164]
            zrangePerTemplate["sbdataExtensionData.dat"]=[0.0, 1.9164]
            zrangePerTemplate["Scd.txt"]=[0.0, 0.5008]
            zrangePerTemplate["StarBurst1.txt"]=[0.0, 0.5008]
            zrangePerTemplate["StarBurst2.txt"]=[0.0, 0.5008]
            zrangePerTemplate["StarBurst3.txt"]=[0.0, 0.5008]
            zrangePerTemplate["vvds_reddestdataExtensionData.dat"]=[0.0, 0.6080]
            zrangePerTemplate["zcosmos_red.txt"]=[1.5109, 3.2469] 
            
        if 1:
            #from Amazed outputs, chi2nc 20160420: ExtendedGalaxyEL2/, for lambda range = [3800, 9459]
            zrangePerTemplate["COMBINE-ave-BX-highblue-AND-Scd.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-BX-highblue-AND-StarBurst1.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-BX-highblue-AND-StarBurst3.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-Lya-abs-AND-Scd.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-Lya-abs-AND-StarBurst1.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-Lya-abs-AND-StarBurst3.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-Lya-emstr-AND-Scd.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-Lya-emstr-AND-StarBurst1.txt"]=[0.0, 6.6]
            zrangePerTemplate["COMBINE-ave-Lya-emstr-AND-StarBurst3.txt"]=[0.0, 6.6]
            zrangePerTemplate["NEW_Im_extended_blue.dat"]=[0.0, 6.6]
            zrangePerTemplate["NEW_Im_extended.dat"]=[0.0, 6.6]
            zrangePerTemplate["NEW_Sbc_extended.dat"]=[0.0, 6.6]
            
            zrangePerTemplate["E_reddataExtensionData.dat"]=[0.0, 2.082]
            zrangePerTemplate["Rebinned_EW-SB2extended.txt"]=[0.0, 3.607]
            zrangePerTemplate["Rebinned_NEW-E-extendeddataExtensionData.txt"]=[0.0, 4.517]
            zrangePerTemplate["vvds_reddestdataExtensionData.dat"]=[0.0, 0.696]
            zrangePerTemplate["COMBINE-zcosmos-red-0Extended-AND-BulgedataExtensionData.txt"]=[0.0, 3.75]
            zrangePerTemplate["COMBINE-zcosmos-red-0Extended-AND-EllipticaldataExtensionData.txt"]=[0.0, 3.75]
            zrangePerTemplate["COMBINE-zcosmos-red-0Extended-AND-s0dataExtensionData.txt"]=[0.0, 3.75]
            zrangePerTemplate["COMBINE-zcosmos-red-0Extended-AND-sadataExtensionData.txt"]=[0.0, 3.75]
            zrangePerTemplate["COMBINE-zcosmos-red-0Extended-AND-sbdataExtensionData.txt"]=[0.0, 3.75]
        
        try:
            zrange = zrangePerTemplate[tplTag]
        except:
            print("Can't find zrange for template tag : {}".format(tplTag))
        
        return zrange
        
    def getBestZcandidateWithinZrangePerTpl(self, chi2Type="raw", extremaType="full", enableZrangeFilter=False):
        zbest = []
        meritbest = []
        tplbest = []
        
        nextrema = 1
        n=self.n
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            if extremaType=="amazed":
                redshifts, merits = self.getZCandidatesFromAmazedChi2Extrema(x, chi2Type, nextrema, enableZrangePerTpl=True)
            else:
                z, merit, tpl = self.getBestZinZrange(x, chi2Type, enableZrangeFilter=enableZrangeFilter)
                
            zbest.append(1e6)
            meritbest.append(1e6)
            tplbest.append("-1")
            
            zbest[x] = z
            meritbest[x] = merit
            tplbest[x] = tpl
                    
            print("zbest = {}".format(zbest[x]))
            
        return zbest, meritbest, tplbest
      
    def getClosestZcandidateZrefIndex(self, indice, chi2Type="raw", extremaType="amazed"):
        print("\n")
        print("Spc {}, getClosestZcandidateZrefIndex".format(indice))
        if extremaType=="amazed":
            redshifts, merits = self.getZCandidatesFromAmazedChi2Extrema(indice, chi2Type)
        else:
            redshifts, merits = self.getZCandidatesFromChi2Extrema(indice, chi2Type)
            
        
        print("redshifts ={}".format(redshifts))
                
        zclosest = redshifts[0]
        indexclosest = 0
        for k in range(len(redshifts)):
            if abs(self.list[indice].zref-redshifts[k]) < abs(self.list[indice].zref-zclosest):
                zclosest = redshifts[k]
                indexclosest = k
        print("iclosest = {}, zref = {}, z = redshifts ={}".format(indexclosest, self.list[indice].zref, zclosest))
                        
        return indexclosest        
            
    def getClosestZcandidatesZrefList(self, chi2Type="raw", extremaType="amazed", enablePlot = False, enableExport = False, nextrema = 5):
        zclosest = []
        zabsdiff = []
        zrank = []
        #n=10
        n=self.n
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            if extremaType=="amazed":
                try:
                    redshifts, merits = self.getZCandidatesFromAmazedChi2Extrema(x, chi2Type, nextrema)
                except Exception as e:
                    print(e)
                    redshifts = []
                    merits = []
                    
            else:
                redshifts, merits = self.getZCandidatesFromChi2Extrema(x, chi2Type, nextrema)
                
            if len(redshifts)<1:
               print("no redshifts found...") 
               print("redshifts ={}".format(redshifts))
            #print("\n\n")
            #print("zref = {}".format(self.list[x].zref))
            zclosest.append(1e6)
            zabsdiff.append(1e6)
            zrank.append(1e6)
            for idx,z in enumerate(redshifts):
                if abs(self.list[x].zref-z) < abs(self.list[x].zref-zclosest[x]):
                    zclosest[x] = z
                    zabsdiff[x] = abs(self.list[x].zref-z)/(1.0+self.list[x].zref)
                    zrank[x] = idx
            
            #print("zabsdiff = {}".format(zabsdiff[x]))
            #plt.plot(zabsdiff)
            #plt.show()
        #print zabsdiff
            
        #prepare subfolder
        if enableExport:            
            outdir = os.path.join(self.analysisoutputdir, "zclosest")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
              
        #export raw data 
        if enableExport:        
            outFile = os.path.join(outdir, 'zclosest_{}_nextrema{}.txt'.format(chi2Type, nextrema))
            outFileFailures = os.path.join(outdir, 'zclosest_{}_nextrema{}_failures5em3.txt'.format(chi2Type, nextrema))
            diffthresFailures = 5e-3
            if not outFile == "":
                # convert to percentage, and save in file 
                fout = open( outFile, "w" ) 
                foutFailures = open( outFileFailures, "w" ) 
                
                fout.write( "idx" + '\t' + "spc" + '\t' + 'zclosest' + '\t' + "rank" + '\t' + "zrelerr"  + '\n')
                foutFailures.write( "idx" + '\t' + "spc" + '\t' + 'zclosest' + '\t' + "rank" + '\t' + "zrelerr"  + '\n')
                for ex in range(0, n):
                    outStr = str(ex) + '\t' + str(self.list[ex].name) + '\t' +str(zclosest[ex]) + '\t' +str(zrank[ex]) + '\t' +str(zabsdiff[ex])
                    #print outStr
                    fout.write( outStr  + '\n')
                    if zabsdiff[ex]>=diffthresFailures:
                        foutFailures.write( outStr  + '\n')
                fout.close()
                foutFailures.close()
            if not outFile == "":
                # convert to percentage, and save in file 
                fout = open( outFile, "w" ) 
                fout.write( "idx" + '\t' + "spc" + '\t' + 'zclosest' + '\t' + "rank" + '\t' + "abszdiff"  + '\n')
                for ex in range(0, n):
                    outStr = str(ex) + '\t' + str(self.list[ex].name) + '\t' +str(zclosest[ex]) + '\t' +str(zrank[ex]) + '\t' +str(zabsdiff[ex])
                    #print outStr
                    fout.write( outStr  + '\n')
                fout.close()
                    
        #export hist data   
        if enablePlot or enableExport:
            fig = plt.figure('closest_zcandidates_found_n{}_chi2{}_extrema{}_histogram'.format(len(zabsdiff), chi2Type, nextrema))
            ax = fig.add_subplot(111)
            #ax.plot(vectErrorBins, yVectErrorBins)
            if 0:            
                plt.hist(zabsdiff, bins=500, normed=1, histtype='stepfilled', cumulative=True)
            else:
                vectErrorBins = np.logspace(-5, 1, 50, endpoint=True)
                ybins = lstats.exportHistogram(yvect=zabsdiff, bins=vectErrorBins, outFile="")
                plt.bar(vectErrorBins, ybins)
            #plt.semilogx()
            plt.xlim([1e-5, 10])
            plt.ylim([0, 100])
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            plt.ylabel('Cumulative Histogram')
            plt.xlabel('abs(zref - closest z candidate) / (1+zref)')
            name1 = "Closest Z candidates found in chi2{}\nN={}, diffthres={}, N extrema ={}".format(chi2Type,len(zabsdiff), self.diffthreshold, nextrema)
            plt.title(name1)
            
            if enableExport:
                outFigFile = os.path.join(outdir, 'closestz_{}_nextrema{}_hist.png'.format(chi2Type, nextrema))
                plt.savefig( outFigFile, bbox_inches='tight')
            if enablePlot:
                plt.show()  


        if enableExport:        
            # ******* large bins histogram
            vectErrorBins = [0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 1.0, 10.0]
            print 'the rough bins are: ' + str(vectErrorBins)
            foutpath = os.path.join(outdir, 'closestz_{}_nextrema{}_stats_brief.txt'.format(chi2Type, nextrema))
            lstats.exportHistogram(zabsdiff, vectErrorBins, foutpath)
        
            # ******* fine bins histogram
            vectErrorBins = np.logspace(-5, 1, 50, endpoint=True)
            print 'the fine bins are: ' + str(vectErrorBins)
            foutpath = os.path.join(outdir, 'closestz_{}_nextrema{}_stats.txt'.format(chi2Type, nextrema))
            lstats.exportHistogram(zabsdiff, vectErrorBins, foutpath)
                                            
        return zabsdiff
        
    def getClosestZcandidateRelativePositionList(self, chi2Type="raw", extremaType="amazed", enablePlot=True, enableExport=False, nextrema = 5):
        zclosest = []
        mclosest = []
        zabsdiff = []
        zrelativemerit = []
        zrelativeredshift = []
        n = self.n
        #n = 10#self.n
        
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            if extremaType=="amazed":
                redshifts, merits = self.getZCandidatesFromAmazedChi2Extrema(x, chi2Type, nextrema)
            else:
                redshifts, merits = self.getZCandidatesFromChi2Extrema(x, chi2Type, nextrema)
                
            #print("redshifts ={}".format(redshifts))
            #print("\n\n")
            #print("zref = {}".format(self.list[x].zref))
            zclosest.append(1e6)
            mclosest.append(1e6)
            if chi2Type=="corr":
                mclosest[x] = -1e6
            zabsdiff.append(1e6)            
            thres = 0.005          
            #thres = 0.01 #for continuum
            for k in range(len(redshifts)):      
                z = redshifts[k]
                conditionHighestMerit = merits[k] < mclosest[x];
                if chi2Type=="corr":
                    conditionHighestMerit = merits[k] > mclosest[x];
                if abs(self.list[x].zref- redshifts[k]) < thres and conditionHighestMerit:
                    zclosest[x] = z
                    mclosest[x] = merits[k]
                    zabsdiff[x] = abs(self.list[x].zref-z)
            #print("zabsdiff = {}".format(zabsdiff))
                    
            z2 = []
            m2 = []
            for k in range(len(redshifts)):                
                if not abs(self.list[x].zref- redshifts[k]) < thres:
                    z2.append(redshifts[k])
                    m2.append(merits[k])
                    
            #print("m2 = {}".format(m2))
            if len(m2) >0:
                sortId=np.argsort(m2)
                #redshifts = [redshiftsUnsorted[b] for b in sortId]
                redshiftssorted = [z2[b] for b in sortId]
                meritssorted = [m2[b] for b in sortId]
                if chi2Type=="corr":
                    zrelativeredshift.append(zclosest[x] - redshiftssorted[len(redshiftssorted)-1])
                    zrelativemerit.append(mclosest[x] - meritssorted[len(meritssorted)-1])
                else:
                    zrelativeredshift.append(zclosest[x] - redshiftssorted[0])
                    zrelativemerit.append(mclosest[x] - meritssorted[0])
            else:
                if chi2Type=="corr":
                    zrelativeredshift.append(1e6)
                    zrelativemerit.append(1e6)
                else:
                    zrelativeredshift.append(1e6)
                    zrelativemerit.append(1e6)
            print("zrelativeredshift = {}".format(zrelativeredshift))
            print("zrelativemerit = {}".format(zrelativemerit))
            #print("zabsdiff = {}".format(zabsdiff[x]))
            #plt.plot(zabsdiff)
            #plt.show()
        #print zabsdiff
        if enableExport:            
            outdir = os.path.join(self.analysisoutputdir, "zrelpos")
            if not os.path.exists(outdir):
                os.makedirs(outdir)  
        
        if enableExport or enablePlot:
            fig = plt.figure('relative_merit_closest_zcandidates_found_chi2{}'.format(chi2Type))
            ax = fig.add_subplot(111)            
            for x in range(0,n):
                if zrelativemerit[x]<1e5:
                    ax.plot(x, zrelativemerit[x], 'bx')
            #plt.hist(zrelativemerit, 500, normed=1, histtype='stepfilled', cumulative=True)
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            #ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            #plt.ylabel('Cumulative Histogram')
            plt.xlabel('merit(closest) - merit(next best extrema)')
            name1 = "merit(closest) - merit(next best extrema) in chi2{} (extrema)".format(chi2Type)
            plt.title(name1)
            if enableExport:
                outFile = os.path.join(outdir, 'zrelpos_{}_nextrema{}.txt'.format(chi2Type, nextrema))
                if not outFile == "":
                    # convert to percentage, and save in file 
                    fout = open( outFile, "w" )  
                    for ex in range(0, n):
                        outStr = str(ex) + '\t' + str(self.list[ex].name) + '\t' +str(zrelativeredshift[ex]) + '\t' +str(zrelativemerit[ex])
                        #print outStr
                        fout.write( outStr  + '\n')
                    fout.close()       
                
                outFigFile = os.path.join(outdir, 'zrelpos_{}_nextrema{}.png'.format(chi2Type, nextrema))
                plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
            if enablePlot:
                plt.show()
            
            
        if enablePlot or enableExport:
            fig = plt.figure('relative_merit_closest_zcandidates_found_chi2{}_histogram'.format(chi2Type))
            ax = fig.add_subplot(111)
            #ax.plot(zrelativemerit)
            nbins = 50
            vectErrorBins = np.linspace(-99, 100, nbins, endpoint=True)
            #ax.plot(zrelativemerit)
            #mybins0=[-2, -1, -0.8, -0.6, -0.4, -0.2 , 0.0, 0.2, 0.4, 0.6, 0.8, 1, 2]
            #mybins=[a*10 for a in mybins0]
            mybins = vectErrorBins
            #mybins0=[-2, -1, -0.8, -0.6, -0.4, -0.2 , 0.0, 0.2, 0.4, 0.6, 0.8, 1, 2]
            #mybins=[a*100 for a in mybins0]
            #print("mybins={}".format(mybins))
            centerbins = [(mybins[k]+mybins[k+1])/2.0 for k in range(len(mybins)-1)]
            widthbins = [(mybins[k+1]-mybins[k]) for k in range(len(mybins)-1)]
            #plt.hist(zrelativemerit, bins=mybins, normed=0, histtype='stepfilled', cumulative=False)
            ybins, bin_edges = np.histogram(zrelativemerit, bins=mybins)
            #ybins = lstats.exportHistogram(yvect=zrelativemerit, bins=mybins, outFile="")
            nbins = len(ybins)
            #print len(ybins)
            #plt.bar(mybins[:-1], ybins, widthbins)#, width = 1)
            plt.plot(mybins[:-1], ybins)
            plt.xlim(min(mybins), max(mybins))
            #plt.semilogx()
            #plt.xlim([-2, 2])
            #plt.ylim([0, 100])            
            
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            #ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            plt.ylabel('Histogram')
            plt.xlabel('merit(closest) - merit(best extrema)')
            name1 = "merit(closest) - merit(best extrema) in chi2{} (extrema)".format(chi2Type)
            plt.title(name1)
            if enableExport:
                outFigFile = os.path.join(outdir, 'zrelpos_{}_nextrema{}_hist.png'.format(chi2Type, nextrema))
                plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
            if enablePlot:
                plt.show() 
                 
            
            if enableExport:
                outFile = os.path.join(outdir, 'zrelpos_{}_nextrema{}_hist.txt'.format(chi2Type, nextrema))
                if not outFile == "":
                    # convert to percentage, and save in file 
                    fout = open( outFile, "w" )  
                    for ex in range(0, nbins):
                        outStr = str(ex) + '\t' + str(centerbins[ex]) + '\t' +str(ybins[ex])
                        print outStr
                        fout.write( outStr  + '\n')
                    fout.close()
                    
            print("\n\nzrelpos_{}_nextrema{}_hist.txt".format(chi2Type, nextrema))
            for ex in range(0, nbins):
                outStr = str(ex) + '\t' + str(centerbins[ex]) + '\t' +str(ybins[ex])
                print outStr                           
        return zrelativemerit
        
    def getSecondExtremaZcandidateRelativePositionList(self, chi2Type="raw", extremaType="amazed", enablePlot=True, enableExport=False, nextrema=10):
        zrelativemerit = []
        diffvect = []
        n = self.n
        #n = 10#self.n
                
        
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            if extremaType=="amazed":
                redshifts, merits = self.getZCandidatesFromAmazedChi2Extrema(x, chi2Type, nextrema)
            else:
                redshifts, merits = self.getZCandidatesFromChi2Extrema(x, chi2Type, nextrema)

            print("redshifts: {}".format(redshifts))       
            print("merits: {}".format(merits))       
            
            z2 = []
            m2 = []
            z2.append(redshifts[0])
            m2.append(merits[0])
            thres = 0.005
            for k in range(1, len(redshifts)):                
                if not abs(redshifts[0] - redshifts[k]) < thres:
                    z2.append(redshifts[k])
                    m2.append(merits[k])
                   
            #print("z2: {}".format(z2))       
            #print("m2: {}".format(m2))
            
            #sortId=np.argsort(m2)
            #meritssorted = [m2[b] for b in sortId]
            zrelativemerit.append(m2[1]-m2[0])
            diffvect.append(self.list[x].zdiff)
        
        if enableExport:            
            outdir = os.path.join(self.analysisoutputdir, "zrelpos_2stExtrema")
            if not os.path.exists(outdir):
                os.makedirs(outdir)   
     
        
        if enablePlot:
            fig = plt.figure('relative_merit_2stExtrema_zcandidates_found_chi2{}'.format(chi2Type))
            ax = fig.add_subplot(111)            
            for x in range(0,n):
                if zrelativemerit[x]<1e5:
                    ax.plot(x, zrelativemerit[x], 'bx')
            #plt.hist(zrelativemerit, 500, normed=1, histtype='stepfilled', cumulative=True)
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            #ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            plt.ylabel('merit(2stExtrema) - merit(best extrema)')
            plt.xlabel('spc id')
            name1 = "merit(2stExtrema) - merit(best extrema) in chi2{} (extrema)".format(chi2Type)
            plt.title(name1)
            if enableExport: 
                outFigFile = os.path.join(outdir, 'relative_merit_2stExtrema_zcandidates_found_chi2{}_nextrema{}_hist.png'.format(chi2Type, nextrema))
                plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
            plt.show()
            
            
        if enablePlot or enableExport:           
            fig = plt.figure('relative_merit_2stExtrema_zcandidates_found_chi2{}_histogram'.format(chi2Type))
            ax = fig.add_subplot(111)
            nbins = 50
            #vectErrorBins = np.linspace(-100, 99, nbins, endpoint=True)
            #ax.plot(zrelativemerit)
            #mybins0=[-2, -1, -0.8, -0.6, -0.4, -0.2 , 0.0, 0.2, 0.4, 0.6, 0.8, 1, 2]
            #mybins=[a*10 for a in mybins0]
            #mybins = vectErrorBins
            #vectErrorBins = np.logspace(-5, 5, nbins, endpoint=True)
            #mybins=[-1.0*a for a in vectErrorBins]
            vectBins = np.logspace(-2, 4, nbins, endpoint=True)
            mybins = vectBins
            
            
            if enableExport:
                foutpath = os.path.join(outdir, 'zrelpos_merit_2stExtrema_zcandidates_found_chi2{}_histogram.txt'.format(chi2Type))
            else:
                foutpath = ""
            ybins = lstats.exportHistogram(yvect=zrelativemerit, bins=mybins, outFile=foutpath, htype='simple', ytype="count")
            #plt.plot(vectErrorBins, ybins, 1.0/nbins)
            plt.plot(mybins, ybins, '-kx')
            
            #plt.semilogx()
            #plt.xlim([0, 1])
            #plt.ylim([0, 100])
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            plt.ylabel('Count')
            plt.xlabel('merit(2stExtrema) - merit(best extrema)')
            #name1 = "Closest Z candidates found in chi2{}\nN={}, diffthres={}, N extrema ={}".format(chi2Type,len(zabsdiff), self.diffthreshold, nextrema)
            name1 = "Histogram\nmerit(2ndExtrema) - merit(1stExtrema) in chi2{} (extrema)".format(chi2Type)
            plt.title(name1)
            if enableExport:
                outRawTxtFile = os.path.join(outdir, 'zrelpos_{}_2ndextrema.txt'.format(chi2Type))
                fout = open( outRawTxtFile, "w" )  
                for ex in range(0, len(zrelativemerit)):
                    outStr = str(ex) + '\t' + str(self.list[ex].name) + '\t' +str(zrelativemerit[ex])
                    fout.write( outStr  + '\n')
                fout.close()        
                
                outFigFile = os.path.join(outdir, 'zrelpos_{}_nextrema{}_hist.png'.format(chi2Type, nextrema))
                plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
            if enablePlot:
                plt.show() 

                    
            print("\n\nzrelpos_{}_nextrema{}_hist.txt".format(chi2Type, nextrema))
            for ex in range(0, nbins):
                outStr = str(ex) + '\t' + str(mybins[ex]) + '\t' +str(ybins[ex])
                print outStr                           
             
        # ******* plot diff in relativemerit2nd extrema bins
        if 1:
            outFileNoExt = 'stats_versusRelativeMerit2ndExtremum_hist' 
            lstats.PlotAmazedVersusChi2diffRatioHistogram(diffvect, [1.0*a for a in zrelativemerit], outdir, outFileNoExt, 1)
            
            return zrelativemerit

    
    def getChi2CombinationCoeff2DMap(self, indice=0, enablePlot=0, chi2dontloadThres=-1, zthres=0.01, opt_combination=1):       
        spcName = self.list[indice].name
        print("spcname = {}".format(spcName))
        s = rp.ResParser(self.dir)
        
        spc = sp.Spectrum(s.getSpcFullPath(spcName)) 
        nspcSamples = (float(len(spc.xvect)))
        print("nsamples = ".format(nspcSamples))
        
        #chi2type = "raw"
        #z_raw, chi2_raw = s.getChi2MinValList(spcName, chi2type) 
        #max_chi2_all = -1e12
        chi2type = "continuum"
        z_continuum, chi2_continuum = s.getChi2MinValList(spcName, chi2type, dontloadThres=chi2dontloadThres)
        for i in range(len(chi2_continuum)):
            chi2_continuum[i] = chi2_continuum[i]/nspcSamples#/1000.0
#            chi2_continuum[i] = chi2_continuum[i]*nspcSamples
#            if max_chi2_all<chi2_continuum[i] and chi2_continuum[i]<1e20:
#                max_chi2_all = chi2_continuum[i]
        chi2type = "nocontinuum"
        z_nc, chi2_nc = s.getChi2MinValList(spcName, chi2type, dontloadThres=chi2dontloadThres) 
        for i in range(len(chi2_nc)):
            chi2_nc[i] = chi2_nc[i]/nspcSamples
#            chi2_nc[i] = chi2_nc[i]*nspcSamples
#            if max_chi2_all<chi2_nc[i] and chi2_nc[i]<1e20:
#                max_chi2_all = chi2_nc[i]
        chi2type = "linemodel"
        z_lm, chi2_lm = s.getChi2MinValList(spcName, chi2type, dontloadThres=chi2dontloadThres)
        for i in range(len(chi2_lm)):
            chi2_lm[i] = chi2_lm[i]/nspcSamples#/nspcSamples
#            if max_chi2_all<chi2_lm[i] and chi2_lm[i]<1e20:
#                max_chi2_all = chi2_lm[i]
        nz = len(chi2_lm)

        #interpolate continuum on linemodel grid
        f = interpolate.interp1d(z_continuum, chi2_continuum, bounds_error=False, fill_value=1e32, assume_sorted=False)
        chi2_continuum_interp = f(z_lm)
        #interpolate no_continuum on linemodel grid
        f = interpolate.interp1d(z_nc, chi2_nc, bounds_error=False, fill_value=1e32, assume_sorted=False) 
        chi2_nc_interp = f(z_lm)
        
        #converting to numpy arrays
        chi2_lm = np.array(chi2_lm)
        chi2_continuum_interp = np.array(chi2_continuum_interp)
        chi2_nc_interp = np.array(chi2_nc_interp)

        #normalizing by min value
        if 0:
            mini_lm = np.min(chi2_lm) 
            mini_chi2c = np.min(chi2_continuum_interp) 
            mini_chi2nc = np.min(chi2_nc_interp) 
            mini_global = np.min([ mini_lm, mini_chi2c, mini_chi2nc ])
            print("normalizing chisquare curves: mini global = {}".format(mini_global))
            
            chi2_lm = np.divide(chi2_lm, mini_global)
            chi2_continuum_interp = np.divide(chi2_continuum_interp, mini_global)
            chi2_nc_interp = np.divide(chi2_nc_interp, mini_global)
                
        #print max_chi2_all
        
#        for i in range(len(chi2_continuum)):
##            if chi2_continuum[i]>1e12:
##                chi2_continuum[i] *= -1.0
#            print chi2_continuum[i]
#            print max_chi2_all
#            print "(inside exp) = {}".format(-(chi2_continuum[i]/max_chi2_all)/2.0)
#            chi2_continuum[i] = np.exp(-(chi2_continuum[i]/max_chi2_all)/2.0)
#            
#            print "(chi2_continuum[i]) = {}".format(chi2_continuum[i])
#        for i in range(len(chi2_nc)):
##            if chi2_nc[i]>1e12:
##                chi2_nc[i] *= -1.0
#            chi2_nc[i] = np.exp(-(chi2_nc[i]/max_chi2_all)/2.0)
#        for i in range(len(chi2_lm)):
##            if chi2_lm[i]>1e12:
##                chi2_lm[i] *= -1.0
#            #chi2_lm[i] = chi2_lm[i]# /nspcSamples
#            chi2_lm[i] = np.exp(-(chi2_lm[i]/max_chi2_all)/2.0)        
        
        if enablePlot==1:
            #plot chi2
            fig = plt.figure('dtreeb - chi2s')
            ax = fig.add_subplot(111)
            ax.plot(z_lm, chi2_nc_interp, label='no continuum')
            ax.plot(z_lm, chi2_lm, label='linemodel', linestyle = "dashed")
            ax.plot(z_lm, chi2_continuum_interp, label='continuum')
            plt.grid(True)
            plt.legend()
            plt.ylabel('Chi2')
            plt.xlabel('z')
            name1 = "Chi2 for {}".format(spcName)
            plt.title(name1)
            
            
            plt.show()
            
            if 0:
                fig = plt.figure('combined merit')
                ax = fig.add_subplot(111)
                ax.plot(z_nc, merit)
                plt.grid(True) 
                plt.show()
        
        #range for coeffs
        coeff_nc = 75
        coeff_continuum_min = 0.0
        coeff_continuum_max = 99.0
        coeff_continuum_step = 2.0
        n_continuum = int((coeff_continuum_max-coeff_continuum_min+1)/float(coeff_continuum_step))
        
        coeff_lm_min = 0.0
        coeff_lm_max = 99.0
        coeff_lm_step = 2.0
        n_lm = int((coeff_lm_max-coeff_lm_min+1)/float(coeff_lm_step))
        merit = np.zeros((nz))
        _map = np.zeros((n_continuum, n_lm))
        
        #print("n_continuum = {}".format(n_continuum))
        #print("n_lm = {}".format(n_lm))
        opt_combine = opt_combination #0=lincomb, 1=bayescomb
        for icontinuum in range(n_continuum):
            print "."
            coeff_continuum = coeff_continuum_min+icontinuum*coeff_continuum_step
            for ilinemodel in range(n_lm):
                coeff_lm = coeff_lm_min+ilinemodel*coeff_lm_step
                #print("this solution uses : coeff_continuum={}, coeff_lm={}".format(coeff_continuum, coeff_lm) )
                if opt_combine==0: #lincomb
                    merit = coeff_nc*chi2_nc_interp + coeff_continuum*chi2_continuum_interp + coeff_lm*chi2_lm
                else:   
                    if 1: #log scale 0-100, goes to 1e-7 to 1.0
                        vectTransform = np.logspace(-10, 0, 101)-10**(-10)
                        prior_nc = vectTransform[int(coeff_nc)]
                        prior_c = vectTransform[int(coeff_continuum)]
                        prior_lm = vectTransform[int(coeff_lm)]  
                    if 0: #log scale
                        logVect = np.logspace(-1.25, 0.025, 101)-10**(-1.25)
                        prior_nc = logVect[int(coeff_nc)]
                        prior_c = logVect[int(coeff_continuum)]/1e5
                        prior_lm = logVect[int(coeff_lm)]
                    if 0: #lin scale
                        prior_nc = coeff_nc/100.0;
                        prior_c = coeff_continuum/100.0; 
                        prior_lm = coeff_lm/100.0;
                    logCoeffLm = -2*np.log(prior_lm)
                    logCoeffC = -2*np.log(prior_c)
                    logCoeffNC = -2*np.log(prior_nc)
                        
                    mini_lm = np.min(chi2_lm+logCoeffLm) 
                    mini_chi2c = np.min(chi2_continuum_interp+logCoeffC) 
                    mini_chi2nc = np.min(chi2_nc_interp+logCoeffNC) 
                    mini_global = np.min([ mini_lm, mini_chi2c, mini_chi2nc ])
                    #print("normalizing chisquare curves: mini global = {}".format(mini_global))
                    
                    valDenom = chi2_lm + logCoeffLm - mini_global
                    likelihood_lm = np.exp(-valDenom/2.0)
                    valDenom = chi2_nc_interp + logCoeffNC - mini_global
                    likelihood_chi2nc= np.exp(-valDenom/2.0)
                    valDenom = chi2_continuum_interp + logCoeffC - mini_global
                    likelihood_chi2c = np.exp(-valDenom/2.0)
                    merit = -( likelihood_lm + likelihood_chi2nc + likelihood_chi2c )

                    
                izbest = np.argmin(merit)
                zbest = z_nc[izbest]
                if(np.abs(zbest - self.list[indice].zref ) < zthres):
                    _map[icontinuum,ilinemodel] = 1
                    print("this solution works ilm={}, ic={}  : z={}, zref={}".format(ilinemodel, icontinuum, zbest, self.list[indice].zref) )
                else:          
                    _map[icontinuum,ilinemodel] = 0
                    print("this solution fails ilm={}, ic={}  : z={}, zref={}".format(ilinemodel, icontinuum, zbest, self.list[indice].zref) )
                    
        return _map
        
    
    def getBayesCombinationCoeff(self, indice=0, enablePlot=0, chi2dontloadThres=-1, zthres=0.01):       
        spcName = self.list[indice].name
        print("spcname = {}".format(spcName))
        respa = rp.ResParser(self.dir)
        

        candid_zvalCandidates, candid_displayParamsBundle = respa.getAutoCandidatesList(spcName)
        lm_candid_displayParamsBundleIdx = 0
        chi2_candid_displayParamsBundleIdx = 1

        #get dtranspose
        filepaths, filenames = respa.getAutoChi2FullPath(spcName)
        filepath = filepaths[1] #1=linemodel
        print("Found chi2 filepath: {}".format(filepath))
        if not os.path.exists(filepath):
            print("Problem while retrieving chi2 filepath.. using: {}".format(filepath))
            stoooooop
        else:
            print("using Chi2 file path : ".format(filepath))
            chi2 = chisq.ResultChisquare(filepath)
            dtransposeD = chi2.amazed_dTransposeD[0]
            dtransposeDNoContinuum = chi2.amazed_dTransposeDNoContinuum[0] #chi2.getFluxMedian()        
        
        #converting to numpy arrays
        z_lm = np.array(candid_zvalCandidates)
        nz = len(z_lm)
        print("z_lm = {}".format(z_lm))
        chi2_lm = np.array(candid_displayParamsBundle[lm_candid_displayParamsBundleIdx]['merits'])
        print("chi2_lm = {}".format(chi2_lm))
        chi2_continuum = np.array(candid_displayParamsBundle[chi2_candid_displayParamsBundleIdx]['merits'])
        print("chi2_continuum = {}".format(chi2_continuum))       
        
        if enablePlot==1:
            #plot chi2
            fig = plt.figure('dtreeb - chi2s')
            ax = fig.add_subplot(111)
            ax.plot(z_lm, chi2_lm, label='linemodel', linestyle = "dashed")
            ax.plot(z_lm, chi2_continuum, label='chi2')
            plt.grid(True)
            plt.legend()
            plt.ylabel('Chi2')
            plt.xlabel('z')
            name1 = "Chi2 for {}".format(spcName)
            plt.title(name1)
            plt.show()
         
        #range for coeffs
        #coeff for sigma compensation (=chi2cSigmaCoeff)
        #coeff_sigmacoeff_min = 0.5
        #coeff_sigmacoeff_max = 1.0
        #coeff_sigmacoeff_step = .25
        #n_sigmacoeff = int((coeff_sigmacoeff_max-coeff_sigmacoeff_min+coeff_sigmacoeff_step)/float(coeff_sigmacoeff_step))
        #coeff_sigmacoeff = np.linspace(0.5, 2.5, 5) 
        coeff_sigmacoeff = [0.1, 0.5, 1.0, 5.0, 10., 100., 500.] 
        n_sigmacoeff = len(coeff_sigmacoeff)
        print("n_sigmacoeff = {}".format(n_sigmacoeff))
        print("coeff_sigmacoeff = {}".format(coeff_sigmacoeff))
        
        #coeff for bayes combination chi2cCoeff
        #coeff_bcoeff_min = -1e4
        #coeff_bcoeff_max = 0.0
        #coeff_bcoeff_step = 100.0
        #n_bcoeff= int((coeff_bcoeff_max-coeff_bcoeff_min+coeff_bcoeff_step)/float(abs(coeff_bcoeff_step)))
        #coeff_bcoeff = 1e4 - (1e4-np.logspace(1.0, 4.0, 20))
        #coeff_bcoeff = np.logspace(-4.0, 0.0, 20)
        coeff_bcoeff = []
        neg_list = - 6.*(np.logspace(1, 3, 10))
        for a in neg_list[::-1]:
            coeff_bcoeff.append(a)
        coeff_bcoeff.append(0.0)
        pos_list = 100.*np.logspace(1, 3, 30)
        for a in pos_list:
            coeff_bcoeff.append(a)
        
        #coeff_bcoeff = np.linspace(0.0, 1e4, 50)
        n_bcoeff = len(coeff_bcoeff)
        print("n_bcoeff = {}".format(n_bcoeff))
        print("coeff_bcoeff = {}".format(coeff_bcoeff))
        
        _map = np.zeros((n_sigmacoeff, n_bcoeff))
        
        merit = np.zeros((nz))

        
        mini_lm_raw = np.min(chi2_lm) 
        mini_chi2_raw = np.min(chi2_continuum) 

        #loop on the sigma and b coeffs to estimate zcalc for each of these coeffs combination
        verbose = 0
        for i_s, coeff_sigma in enumerate(coeff_sigmacoeff):
            print "."
            for i_b, coeff_b in enumerate(coeff_bcoeff):
                coeff_b_negative = -coeff_b
                #coeff_b_negative = -dtransposeD*coeff_b
                if verbose:
                    print("this solution uses : coeff_sigma={}, coeff_b={:.1f}".format(coeff_sigma, coeff_b) )
                #continue
                
                if 1: 
                    mini_lm = np.min(chi2_lm) 
                    mini_chi2 = np.min(chi2_continuum*coeff_sigma+coeff_b_negative) 
                    mini_global = np.min([ mini_lm, mini_chi2 ])
                    if verbose:
                        print("normalizing chisquare curves: mini global = {}".format(mini_global))
                    
                    valDenom = chi2_lm - mini_global
                    likelihood_lm = np.exp(-valDenom/2.0)
                    valDenom = chi2_continuum + coeff_b_negative - mini_global
                    likelihood_chi2 = np.exp(-valDenom/2.0)
                    merit = -np.log( likelihood_lm + likelihood_chi2 )
                    if verbose:
                        print("combined merit = {}".format(merit))
                                        
                izbest = np.argmin(merit)
                if verbose:
                    print("izbest = {}".format(izbest))
                zbest = z_lm[izbest]
                if(np.abs(zbest - self.list[indice].zref )/(1.+self.list[indice].zref) < zthres):
                    _map[i_s,i_b] = 1
                    if verbose:
                        print("this solution works i_s={}, i_b={}  : z={}, zref={}".format(i_s, i_b, zbest, self.list[indice].zref) )
                else:          
                    _map[i_s,i_b] = 0
                    if verbose:
                        print("this solution fails i_s={}, i_b={}  : z={}, zref={}".format(i_s, i_b, zbest, self.list[indice].zref) )
                
                    
        return _map, coeff_sigmacoeff, coeff_bcoeff, dtransposeDNoContinuum, dtransposeD
    
    def plotReducedZcandidates(self, chi2Type="raw", extremaType="amazed"):
        zrchi, mrchi = self.getReducedZcandidates(chi2Type="raw")
        zrchinocontinuum, mrchinocontinuum = self.getReducedZcandidates(chi2Type="nocontinuum")
        zrcorr, mrcorr = self.getReducedZcandidates(chi2Type="corr")
        n = 30
        if True:
            fig = plt.figure('reduced_zcandidates_found_chi2{}'.format(chi2Type))
            ax = fig.add_subplot(111)
        
            thres = 0.005
            for x in range(0,n):
                
                for nres in range(3):
                    if nres == 0:
                        mc = mrchi
                        zc = zrchi
                        c = "bx"
                        offset = 0.0
                    elif nres == 1:
                        mc = mrchinocontinuum
                        zc = zrchinocontinuum
                        c = "kx"
                        offset = 0.25
                    elif nres == 2:
                        mc = -np.copy(mrcorr)+1
                        zc = zrcorr
                        c = "gx"
                        offset = 0.5
                        
                    npts = len(mc[x])
                    for k in range(npts):
                        #print("plotting point {}/{}".format(k,npts))
                        if not abs(self.list[x].zref-zc[x][k]) < thres:
                            ax.plot(x + offset, mc[x][k], c)
                        else:
                            ax.plot(x + offset + 0.05, mc[x][k], 'ro')                        
            
            #plt.hist(zabsdiff, 500, normed=1, histtype='stepfilled', cumulative=True)
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            #ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            #plt.ylabel('Cumulative Histogram')
            #plt.xlabel('abs(zref - closest z candidate)')
            #name1 = "Closest Z candidates found in chi2{} (extrema)".format(chi2Type)
            #plt.title(name1)
            plt.show()
    
    def getReducedZcandidates(self, chi2Type="raw", extremaType="amazed", enablePlot = False):
        mc = []
        zz = []
        n = 30#self.n
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            if extremaType=="amazed":
                redshifts, merits = self.getZCandidatesFromAmazedChi2Extrema(x, chi2Type)
            else:
                redshifts, merits = self.getZCandidatesFromChi2Extrema(x, chi2Type)
                
            print("redshifts ={}".format(redshifts))
            print("merits ={}".format(merits))
            
            zz.append(redshifts)
            #print("\n\n")
            #print("zref = {}".format(self.list[x].zref))
            #merits = 
            mmax = np.max(merits)
            mmin = np.min(merits)
            meritsc = []
            for k in range(len(merits)):
                meritsc.append((merits[k]-mmin)/(mmax-mmin))
            #print("zabsdiff = {}".format(zabsdiff[x]))
            #plt.plot(zabsdiff)
            #plt.show()
            mc.append(meritsc)
            
        
        #print zabsdiff
        if enablePlot:
            fig = plt.figure('reduced_zcandidates_found_chi2{}'.format(chi2Type))
            ax = fig.add_subplot(111)
        
            thres = 0.005
            for x in range(0,n):
                #xvect = np.ones(len(zc[x]))*x
                #kbest = self.getClosestZcandidateZrefIndex(x, chi2Type, extremaType)
                npts = len(mc[x])
                for k in range(npts):
                    #print("plotting point {}/{}".format(k,npts))
                    if not abs(self.list[x].zref-zz[x][k]) < thres:
                        ax.plot(x, mc[x][k], 'bx')
                    else:
                        ax.plot(x+0.05, mc[x][k], 'ro')                        
            #plt.hist(zabsdiff, 500, normed=1, histtype='stepfilled', cumulative=True)
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            #ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            #plt.ylabel('Cumulative Histogram')
            #plt.xlabel('abs(zref - closest z candidate)')
            #name1 = "Closest Z candidates found in chi2{} (extrema)".format(chi2Type)
            #plt.title(name1)
            plt.show()                                      
        return zz, mc

    def getPeakDetectionComparison(self, refreslist, enablePlot=True, enableExport=True):
        '''
        This function looks if all the peaks detected in the refresdir matching/
        result are also detected in the current matching result
        '''
        s = rp.ResParser(self.dir)
        sref = rp.ResParser(refreslist.dir)
        #print(s) 
        samedetectionrate = []
        npeaksdetected = []
        
        n = self.n
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            spcName = self.list[x].name
            if spcName.find("-F_")!=-1:
                spcRefName = spcName.replace("-F_", "-TF_")
            elif spcName.find("FILTERED_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version")!=-1:
                spcRefName = spcName.replace("FILTERED_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version", "EZ_fits-W-TF_")
            else:
                print("the line detection lists are not the same. Spectrum Names differ from one resList to another... aborting")
                return
            print("spcRefName: {}".format(spcRefName))
            ctlgpath = s.getDetectedPeakCatalogPath(spcName)
            catalog = dpeakctlg.DetectedPeakCatalog(ctlgpath)
            ctlgrefpath = sref.getDetectedPeakCatalogPath(spcRefName)
            catalogref = dpeakctlg.DetectedPeakCatalog(ctlgrefpath)           
            samedetectionrate.append(0)
            for kref in range(catalogref.n):
                centerindex = catalogref.list[kref].inf + (catalogref.list[kref].sup-catalogref.list[kref].inf)/2.0
                for k in range(catalog.n):
                    if centerindex > catalog.list[k].inf and centerindex < catalog.list[k].sup :
                        samedetectionrate[x]+=1
                        break
            samedetectionrate[x] = float(samedetectionrate[x])/float(catalogref.n)
            npeaksdetected.append(catalog.n)

        nmoy = np.mean(npeaksdetected) 
        nmax = np.max(npeaksdetected) 
        nmin = np.min(npeaksdetected) 
        print("detectedcatalog stats, npeaksmoy = {}, npeaksmax = {}, nlinemin = {}".format(nmoy, nmax, nmin))
        samedetectionratemoy = np.mean(samedetectionrate)
        samedetectionratemin = np.min(samedetectionrate)
        samedetectionratemax = np.max(samedetectionrate)
        samedetectionratefailures = np.argwhere(samedetectionrate<1.0)
        print("detectedcatalog comparison with ref, samdetectionrate = {}, samdetectionratemax = {}, samdetectionratemin = {}".format(samedetectionratemoy, samedetectionratemax, samedetectionratemin))
        print("detectedcatalog comparison with ref, samedetectionratefailures = {}".format(samedetectionratefailures))
        
        if enableExport:            
            outdir = os.path.join(self.analysisoutputdir, "linedetection")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                
        if enablePlot:
            if 0:
                fig = plt.figure('N lines found')
                ax = fig.add_subplot(111)
                ax.plot(npeaksdetected)
                plt.show() 
            
            fig = plt.figure('peakdetection_comparison_histogram')
            ax = fig.add_subplot(111)
            #ax.plot(vectErrorBins, yVectErrorBins)
            nbins = 50
            vectErrorBins = np.linspace(0, 1, nbins, endpoint=True)
            if enableExport:
                foutpath = os.path.join(outdir, 'samepeakdetectionrate_stats.txt')
            else:
                foutpath = ""
            ybins = lstats.exportHistogram(yvect=samedetectionrate, bins=vectErrorBins, outFile=foutpath, htype='simple')
            #plt.plot(vectErrorBins, ybins, 1.0/nbins)
            plt.plot(vectErrorBins, ybins)
            
            #plt.semilogx()
            plt.xlim([0, 1])
            plt.ylim([0, 100])
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            #ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            plt.ylabel('Histogram')
            plt.xlabel('Same detection rate')
            #name1 = "Closest Z candidates found in chi2{}\nN={}, diffthres={}, N extrema ={}".format(chi2Type,len(zabsdiff), self.diffthreshold, nextrema)
            #plt.title(name1)
            plt.show() 

        
    
    def getLineDetectionStatsList(self, enablePlot=True, enableExport=False):
        s = rp.ResParser(self.dir)
        #print(s) 
        nlines = []
        nmoy = 0.0
        nmax = -1e6
        nmin = 1e6
        
        n = self.n
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            spcName = self.list[x].name
            ctlgpath = s.getDetectedLineCatalogPath(spcName)
            catalog = dctlg.Catalog(ctlgpath)
            nlines.append(catalog.n)
            nmoy += nlines[x]
            if nmax < nlines[x]:
                nmax = nlines[x]
            if nmin > nlines[x]:
                nmin = nlines[x]
        nmoy /= float(n)
        nstd = 0.0
        for x in range(0,n):
            nstd += (nmoy - nlines[x])**2
        nstd /= float(n)
        nstd = math.sqrt(nstd)
        
        print("detectedcatalog stats, moy = {}, std = {}".format(nmoy, nstd))
        print("detectedcatalog stats, min = {}, max = {}".format(nmin, nmax))
        
        if enablePlot:
            fig = plt.figure('N lines found')
            ax = fig.add_subplot(111)
            ax.plot(nlines)
            plt.show()
 
    def getLineDetectionComparison(self, refreslist, enablePlot=True, enableExport=True):
        '''
        This function looks if all the lines detected in the refresdir matching/
        result are also detected in the current matching result
        '''
        thres = 1.0
        s = rp.ResParser(self.dir)
        sref = rp.ResParser(refreslist.dir)
        #print(s) 
        samedetectionrate = []
        nlinesdetected = []
        
        n = self.n
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            spcName = self.list[x].name
            if spcName.find("-F_")!=-1:
                spcRefName = spcName.replace("-F_", "-TF_")
            elif spcName.find("FILTERED_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version")!=-1:
                spcRefName = spcName.replace("FILTERED_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version", "EZ_fits-W-TF_")
            else:
                print("the line detection lists are not the same. Spectrum Names differ from one resList to another... aborting")
                return
            print("spcRefName: {}".format(spcRefName))
            ctlgpath = s.getDetectedLineCatalogPath(spcName)
            catalog = dctlg.Catalog(ctlgpath)
            ctlgrefpath = sref.getDetectedLineCatalogPath(spcRefName)
            catalogref = dctlg.Catalog(ctlgrefpath)           
            samedetectionrate.append(0)
            for kref in range(catalogref.n):
                for k in range(catalog.n):
                    if abs(catalogref.linelambda[kref] - catalog.linelambda[k]) < thres :
                        samedetectionrate[x]+=1
                        break
            samedetectionrate[x] = float(samedetectionrate[x])/float(catalogref.n)
            nlinesdetected.append(catalog.n)

        nmoy = np.mean(nlinesdetected) 
        nmax = np.max(nlinesdetected) 
        nmin = np.min(nlinesdetected) 
        print("detectedcatalog stats, nlinesmoy = {}, nlinesmax = {}, nlinemin = {}".format(nmoy, nmax, nmin))
        samedetectionratemoy = np.mean(samedetectionrate)
        samedetectionratemin = np.min(samedetectionrate)
        samedetectionratemax = np.max(samedetectionrate)
        samedetectionratefailures = np.argwhere(samedetectionrate<1.0)
        print("detectedcatalog comparison with ref, samdetectionrate = {}, samdetectionratemax = {}, samdetectionratemin = {}".format(samedetectionratemoy, samedetectionratemax, samedetectionratemin))
        print("detectedcatalog comparison with ref, samedetectionratefailures = {}".format(samedetectionratefailures))
        
        if enableExport:            
            outdir = os.path.join(self.analysisoutputdir, "linedetection")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                
        if enablePlot:
            if 0:
                fig = plt.figure('N lines found')
                ax = fig.add_subplot(111)
                ax.plot(nlinesdetected)
                plt.show() 
            
            fig = plt.figure('detection_comparison_histogram')
            ax = fig.add_subplot(111)
            #ax.plot(vectErrorBins, yVectErrorBins)
            vectErrorBins = np.linspace(0, 1, 50, endpoint=True)
            if enableExport:
                foutpath = os.path.join(outdir, 'samedetectionrate_stats.txt')
            else:
                foutpath = ""
            ybins = lstats.exportHistogram(yvect=samedetectionrate, bins=vectErrorBins, outFile=foutpath, htype='simple')
            #plt.bar(vectErrorBins, ybins, 1.0/50.0)
            plt.plot(vectErrorBins, ybins)
            
            #plt.semilogx()
            plt.xlim([0, 1])
            plt.ylim([0, 100])
            ##bar
            #ind = np.arange(len(OY))
            #pp.plot(xvect, yvect, 'x')
            #ax.set_xscale('log')
            plt.grid(True) # Affiche la grille
            #pp.legend(('cos','sin'), 'upper right', shadow = True)
            plt.ylabel('Histogram')
            plt.xlabel('Same detection rate')
            #name1 = "Closest Z candidates found in chi2{}\nN={}, diffthres={}, N extrema ={}".format(chi2Type,len(zabsdiff), self.diffthreshold, nextrema)
            #plt.title(name1)
            plt.show() 

    
    def getLineMatchingStatsList(self, enablePlot=True, enableExport=False):
        s = rp.ResParser(self.dir)
        #print(s) 
        noklist = []
        okmoy = 0.0
        nfaillist = []
        failmoy = 0.0
        failmax = -1e6
        failmin = 1e6
        
        n = self.n
        #n=32  # rule 1, n=47-32, #rule 2, n=69
        
        for x in range(0,n):
            print("\n")
            print("Spc {}/{}".format(x, self.n))
            spcName = self.list[x].name
            print("SpcName = {}".format(spcName))
            matchresultpath = s.getLineMatchingResultPath(spcName)
            matchres = matchsol.MatchingSolution(matchresultpath)
            nok, nfail = matchres.getNCorrectSolution(self.list[x].zref)            
            
            noklist.append(nok)
            okmoy += nok>0
            nfaillist.append(nfail)
            #failratio = float(nfail)/float(nfail+nok)
            failmoy += nfail
            if failmax < nfail:
                failmax = nfail
            if failmin > nfail:
                failmin = nfail
        failmoy /= float(n)
        okmoy /= float(n)
        
        print("matching sol. stats, okmoy = {}".format(okmoy))
        print("matching sol. stats, failmoy = {}".format(failmoy))
        print("matching sol. stats, failmin = {}, failmax = {}".format(failmin, failmax))
        
        if enablePlot:
            fig = plt.figure('LineMatching Solution, nfail')
            ax = fig.add_subplot(111)
            ax.plot(nfaillist)
            plt.grid()
            plt.show()  
            fig = plt.figure('LineMatching Solution, nok')
            ax = fig.add_subplot(111)
            xvect = range(len(nfaillist))
            ax.plot(xvect, noklist, 'b', label="ok")
            #ax.plot(xvect, nfaillist, 'k', label="fail")
            plt.legend()
            plt.grid()
            plt.show()
            
    def getFailuresComparison(self, refreslist, enablePlot=True, enableExport=True):
        '''
        This function compares the failures intersection and differences
        '''
        #s = rp.ResParser(self.dir)
        #sref = rp.ResParser(refreslist.dir)
        #print(s) 
        failures_both = []
        failures_onlyThis = []
        failures_onlyRef = []
        
        zrelerr_onlyThis = []
        zrelerr_onlyRef = []

        removeStrRef = ""
        removeStrRef = "_interleaved"
        #removeStrRef = "_FILT_MGv0_c"
        #removeStrRef = "_FILT3dec2015"
        #removeStrThis = "_FILT_MGv1_c"
        removeStrThis = "_combined"  
        #removeStrThis = ""
        #removeStr2 = "_FILT"              
        WarningKeyStr = raw_input("\n\nWARNING: Removing strings for comparison : ref, removing={} - this, removing={}\nPress any key to continue...".format(removeStrRef, removeStrThis))
        
        n = self.n
        nref = refreslist.n
        print("\nINFO: searching for failures in THIS dataset...")
        for x in range(0,n):
            #print("\n")
            print("Spc {}/{}".format(x+1, self.n))
            spcName = self.list[x].name
            spcName = spcName.replace(removeStrThis, "")
            relzerr = (self.list[x].zcalc-self.list[x].zref)/(1.0+self.list[x].zref)
            spcFound = False;
            for xref in range(0,nref):
                #print("\n")
                #print("SpcRef {}/{}".format(xref+1, nref)) 
                spcNameRef = refreslist.list[xref].name
                spcNameRef = spcNameRef.replace(removeStrRef, "")
                if spcName.find(spcNameRef)!=-1:
                    failures_both.append(spcName)
                    spcFound = True;
                    break
            if not spcFound:
                failures_onlyThis.append(spcName)
                zrelerr_onlyThis.append(relzerr)
                
                
        print("\nINFO: searching for failures in REF dataset...")
        for xref in range(0,nref):
            #print("\n")
            print("Spc {}/{}".format(x+1, self.n))
            spcNameRef = refreslist.list[xref].name
            spcNameRef = spcNameRef.replace(removeStrRef, "")
            relzerrRef = (refreslist.list[xref].zcalc-refreslist.list[xref].zref)/(1.0+refreslist.list[xref].zref)
            spcRefFound = False;
            for x in range(0,n):
                #print("\n")
                #print("SpcRef {}/{}".format(xref+1, nref)) 
                spcName = self.list[x].name
                spcName = spcName.replace(removeStrThis, "")
                if spcName.find(spcNameRef)!=-1:
                    spcRefFound = True;
                    break
            if not spcRefFound:
                failures_onlyRef.append(spcNameRef)
                zrelerr_onlyRef.append(relzerrRef)
                
        print("\n\n")
        print("########## Failures Comparison: ##########")
        print("N same failures = {}".format(len(failures_both)))
        print("\n")
        print("N total failures this = {}".format(n))
        print("N failures only this = {}".format(len(failures_onlyThis)))
        for x in range(0,len(failures_onlyThis)):
            print("\t{}\t{}".format(failures_onlyThis[x], zrelerr_onlyThis[x]))
        print("\n")
        print("N total failures ref = {}".format(nref))
        print("N failures only ref = {}".format(nref-len(failures_both)))
        for x in range(0,len(failures_onlyRef)):
            print("\t{}\t{}".format(failures_onlyRef[x], zrelerr_onlyRef[x]))
        print("\n")
        
        if enableExport:            
            outdir = os.path.join(self.analysisoutputdir, "failures")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            foutpath = os.path.join(outdir, "failures_comparison_with_{}.txt".format(refreslist.name))
            f = open(foutpath, "w")
            sameStr = "N same failures = {}\n\n".format(len(failures_both))
            f.write(sameStr)
            thisStr = "Failures for THIS resultset = {} :\n".format(self.name)
            thisStr = thisStr + "N total failures THIS resultset = {}\n".format(n)
            thisStr = thisStr + "N failures only THIS resultset = {}\n".format(len(failures_onlyThis))        
            for x in range(0,len(failures_onlyThis)):
                thisStr = thisStr +"\t{}\t{}\n".format(failures_onlyThis[x], zrelerr_onlyThis[x])
            f.write(thisStr+"\n")  
            refStr = "Failures for REF resultset = {} :\n".format(refreslist.name)
            refStr = refStr + "N total failures REF resultset = {}\n".format(nref)
            refStr = refStr + "N failures only REF resultset = {}\n".format(nref-len(failures_both))        
            for x in range(0,len(failures_onlyRef)):
                refStr = refStr +"\t{}\t{}\n".format(failures_onlyRef[x], zrelerr_onlyRef[x])
            f.write(refStr+"\n")  
            f.close()

        return

    
    
   
def plotRelativePosClosestZref(resDir, spcName="", nextrema=15, extremaType = 'raw'):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1, opt='brief', spcName=spcName)
    if resList.n <1:
        print('No results loaded...')
        return
    zrelativemeritList = []
    
#    relposraw = resList.getClosestZcandidateRelativePositionList("raw", "amazed", enablePlot=True, enableExport=False, nextrema = nextrema)
#    zrelativemeritList.append(relposraw)
    relposnocont = resList.getClosestZcandidateRelativePositionList(extremaType, "amazed", enablePlot=True, enableExport=True, nextrema = nextrema)
#    zrelativemeritList.append(relposnocont)
#    relposcorr = resList.getClosestZcandidateRelativePositionList("corr", "amazed", enablePlot=True, enableExport=False, nextrema = nextrema)
#    zrelativemeritList.append(relposcorr)
#    relposlinemodel = resList.getClosestZcandidateRelativePositionList("linemodel", "amazed", enablePlot=True, enableExport=True, nextrema = nextrema)
    
#    zrelativemeritList.append(relposnocont)
#    zrelativemeritList.append(relposlinemodel)
    
    n = resList.n
    #n = 10
    if 0:
        fig = plt.figure('relative_merit_closest_zcandidates')
        ax = fig.add_subplot(111) 
        
        cmax = -1e6
        cmin = 1e6
        for x in range(0,n):
            for nres in range(3): 
                if zrelativemeritList[nres][x]<cmin and zrelativemeritList[nres][x]<1e5:
                    cmin = zrelativemeritList[nres][x]
                if zrelativemeritList[nres][x]>cmax and zrelativemeritList[nres][x]<1e5:
                    cmax = zrelativemeritList[nres][x]
                    
        for nres in range(3): 
            if nres == 0:
                c = "bo"
                leg = 'chi2'
            elif nres == 1:
                c = "k^"
                leg = 'chi2nocontinuum'
            elif nres == 2:
                c = "gx"
                leg = 'corr'
            elif nres == 3:
                c = "rx"
                leg = 'linemodel'
            ax.plot(zrelativemeritList[nres], c, label=leg)
                #plt.hist(zrelativemerit, 500, normed=1, histtype='stepfilled', cumulative=True)
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        #ax.set_xscale('log')
        
        plt.ylim([cmin , cmax])
        plt.grid(True) # Affiche la grille
        plt.legend()
        #plt.legend(('raw','no continuum', 'corr'), 'upper right', shadow = True)
        #plt.ylabel('Cumulative Histogram')
        plt.xlabel('merit(closest) - merit(next extrema)')
        #name1 = "merit(closest) - merit(next extrema) in chi2{} (extrema)".format(chi2Type)
        #plt.title(name1)
        plt.show()
        
def plotRelativePosSecondBestExtrema(resDir, diffthres, spcName=""):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName)
    if resList.n <1:
        print('No results loaded...')
        return
    zrelativemeritList = []
    
    nextrema = 10
#    relposraw = resList.getClosestZcandidateRelativePositionList("raw", "amazed", enablePlot=True, enableExport=False, nextrema = nextrema)
#    zrelativemeritList.append(relposraw)
    relposnocont = resList.getClosestZcandidateRelativePositionList("nocontinuum", "amazed", enablePlot=True, enableExport=True, nextrema = nextrema)
    zrelativemeritList.append(relposnocont)
#    relposcorr = resList.getClosestZcandidateRelativePositionList("corr", "amazed", enablePlot=True, enableExport=False, nextrema = nextrema)
#    zrelativemeritList.append(relposcorr)
#    relposlinemodel = resList.getSecondExtremaZcandidateRelativePositionList("linemodel", "amazed", enablePlot=True, enableExport=True, nextrema=nextrema)
#    zrelativemeritList.append(relposlinemodel)
    
    n = resList.n
    #n = 10
    if 0:
        fig = plt.figure('relative_merit_closest_zcandidates')
        ax = fig.add_subplot(111) 
        
        cmax = -1e6
        cmin = 1e6
        for x in range(0,n):
            for nres in range(3): 
                if zrelativemeritList[nres][x]<cmin and zrelativemeritList[nres][x]<1e5:
                    cmin = zrelativemeritList[nres][x]
                if zrelativemeritList[nres][x]>cmax and zrelativemeritList[nres][x]<1e5:
                    cmax = zrelativemeritList[nres][x]
                    
        for nres in range(3): 
            if nres == 0:
                c = "bo"
                leg = 'chi2'
            elif nres == 1:
                c = "k^"
                leg = 'chi2nocontinuum'
            elif nres == 2:
                c = "gx"
                leg = 'corr'
            elif nres == 3:
                c = "rx"
                leg = 'linemodel'
            ax.plot(zrelativemeritList[nres], c, label=leg)
                #plt.hist(zrelativemerit, 500, normed=1, histtype='stepfilled', cumulative=True)
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        #ax.set_xscale('log')
        
        plt.ylim([cmin , cmax])
        plt.grid(True) # Affiche la grille
        plt.legend()
        #plt.legend(('raw','no continuum', 'corr'), 'upper right', shadow = True)
        #plt.ylabel('Cumulative Histogram')
        plt.xlabel('merit(closest) - merit(next extrema)')
        #name1 = "merit(closest) - merit(next extrema) in chi2{} (extrema)".format(chi2Type)
        #plt.title(name1)
        plt.show()   

def exportBestRedshiftWithZRangePerTemplate(resDir, diffthres, chi2Type="linemodeltplshape", spcName="", enableZrangeFilter=False):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName)
    if resList.n <1:
        print('No results loaded...')
        return
    
    bestz, bestmerit, bestTpl = resList.getBestZcandidateWithinZrangePerTpl(chi2Type, "full", enableZrangeFilter=enableZrangeFilter)
    
    enableExport = True
    if enableExport:    
            zrangeFilterSwitchStr = "ON"
            if not enableZrangeFilter:
                zrangeFilterSwitchStr = "OFF"
            outdir = os.path.join(resList.analysisoutputdir, "best_redshift_chi2Type{}_zrangeFilter{}".format( chi2Type, zrangeFilterSwitchStr))
            if not os.path.exists(outdir):
                print("creating outputdir {}".format(outdir))
                os.makedirs(outdir) 
                
            #export new redshift.csv file
            newredshift_fname = os.path.join(outdir, "best_redshift_zrangepertpl.csv")
            text_file = open(newredshift_fname, "w")
            data = "#id   Z   MERIT   TplID   Method\n"
            text_file.write("{}".format(data))
            for k in range(resList.n):    
                data = "{}\t{}\t{}\t{}\t{}\n".format(resList.list[k].name, bestz[k], bestmerit[k], bestTpl[k], "method" )
                text_file.write("{}".format(data))
            text_file.close()
            
def printSourcesInZMagSfrBin(resDir, diffthres, zrefmin, zrefmax, magrefmin, magrefmax, sfrrefmin, sfrrefmax):
    print('\n')
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName="", methodName="", zrefmin=zrefmin, zrefmax=zrefmax, magrefmin=magrefmin, magrefmax=magrefmax, sfrrefmin=sfrrefmin, sfrrefmax=sfrrefmax)
    if resList.n <1:
        print('No results loaded...')
        return
    
    for k in range(resList.n):
        relzerr = (resList.list[k].zcalc-resList.list[k].zref)/(1+resList.list[k].zref)
        print("{:<8}{:<50}\t{:3.4f}\t{:3.4f}".format(k, resList.list[k].name, resList.list[k].zref, relzerr))
    

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in xrange(N+1) ]
    # Return colormap object.
    return mcolors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)
    
def colorbar_index(ncolors, cmap):
    cmap = cmap_discretize(cmap, ncolors)
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(range(ncolors))
      
def plotChi2CombinationCoeff2DMap(resDir, diffthres, spcName="", methodName="", enableExport=True):
    print('using amazed results full path: {0}'.format(resDir))
    zrefmin=0.2
    zrefmax=2.5
    magrefmin=0.0
    magrefmax=23.3
    
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName, methodName=methodName, zrefmin=zrefmin, zrefmax=zrefmax, magrefmin=magrefmin, magrefmax=magrefmax)
    if resList.n <1:
        print('No results loaded...')
        return
    nPtsPerAxis = 50
    cumulcoeffmap = np.zeros((nPtsPerAxis,nPtsPerAxis))
    fullcoeffmap = np.ones((resList.n,nPtsPerAxis,nPtsPerAxis))*-1
    
    ## parameters :
    opt_combination = 0 #0=lincomb, 1=bayescomb 
    zthres = 0.01
    enableDirectPlottingDebugMode = 0
    dont_skip_no_spc = 0
    if not enableDirectPlottingDebugMode: 
        plt.ion()
        chi2dontloadThres = -1#1e12 #better for direct plotting, but use -1 to do the optimization map
        enablePlot = False
    else:
        #plt.ion()
        chi2dontloadThres = 1e12 #better for direct plotting, but use -1 to do the optimization map
        enablePlot = True
    
    #sleep optional
    if False:
        sleep_time_seconds = 3600*2.0
        print("Sleeping for {} seconds...".format(sleep_time_seconds))
        time.sleep(sleep_time_seconds)

    
    for k in range(resList.n):
        print("processing result #{}/{}".format(k+1, resList.n))
        try:
            coeffmap = resList.getChi2CombinationCoeff2DMap(k, enablePlot=enablePlot, chi2dontloadThres=chi2dontloadThres,  zthres=zthres, opt_combination=opt_combination)
        except Exception as e:
            print(e)
            if dont_skip_no_spc:
                stop
            else:
                continue

        print("coeffmap = {}".format(coeffmap))
        fullcoeffmap[k, :, :] = np.copy(coeffmap)
        cumulcoeffmap = np.add(cumulcoeffmap, coeffmap)
        print("cumulcoeffmap = {}".format(cumulcoeffmap))
        #print("coeffmap shape = {}".format(coeffmap.shape)) 
        
        #NOTE: in order to deactivate showing this window at each iteration, use Agg (see imports)
        plt.clf()
        plt.close()
  
        if opt_combination==0:
            comb_name = 'lincomb'
        elif opt_combination == 1:
            comb_name = 'bayescomb'
        
        fig = plt.figure('dtreeb {} coeff map'.format(comb_name), figsize=(9, 8))
        ax = fig.add_subplot(111)
        cmap = plt.get_cmap('RdYlGn') 
        ncolors = 20 #min(20,k+2)
        cmap = cmap_discretize(cmap, ncolors) 

        i = ax.matshow(np.transpose(cumulcoeffmap), interpolation='nearest', aspect='equal', cmap=cmap)

        plt.xlabel('continuum coeff')
        plt.ylabel('lm coeff')
        name1 = "combination coeff map \nzThreshold = {}\n(processed idx={}/{})".format(zthres, k+1, resList.n)
        plt.title(name1)
        #plt.legend(legendz)
        #plt.grid()
        #i = ax.imshow(image, interpolation='nearest')
        #if k==0:
        #    fig.colorbar(i)
        cbar = fig.colorbar(i)
        cbmax = np.amax(cumulcoeffmap)
        i.set_clim(vmin=cbmax-ncolors-0.5, vmax=cbmax+0.5)
        cbar.set_ticks(np.linspace(cbmax-ncolors, cbmax, ncolors))
        cbar.set_ticklabels(range(int(cbmax-ncolors), int(cbmax)))

        #colorbar_index(ncolors=cbmax, cmap=cmap)   
        
        #plt.clim(0,k+1)
        #plt.clim(cmin,cmin*2.0)
        plt.draw()
        
       
        if enableExport:    
            tag = "dtreeb_{}_coeffmap_zthres{}_method{}".format(comb_name, zthres, methodName)
            outdir = os.path.join(resList.analysisoutputdir, tag)
            if not os.path.exists(outdir):
                print("creating outputdir {}".format(outdir))
                os.makedirs(outdir)   
            outFigFile = os.path.join(outdir, 'cumulcoeffmap.png')
            plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
            outtxtFile = os.path.join(outdir, 'cumulcoeffmap_{}.txt'.format(tag))
            np.savetxt(outtxtFile, cumulcoeffmap)
            outFullMapFile = os.path.join(outdir, 'fullcoeffmap_{}.dat'.format(tag))
            np.save(outFullMapFile, fullcoeffmap)

def estimateCombinationCoeffMap(resDir, diffthres, spcName="", enableExport=True):
    #options
    dont_skip_no_spc = 1
    enable_intermediate_plots = 0
    zthres=0.001
    
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName)
    if resList.n <1:
        print('No results loaded...')
        return   
        
 
    success_begin_points = []
    success_end_points = []
    success_range = [] 
    coeffMapCumulative = None
    #plt.ion()
    #nres = 40
    nres = resList.n
    for k in range(nres):
        print("processing result #{}/{}".format(k+1, resList.n))
        print("zref ={}".format(resList.list[k].zref))
        try:
            #coeffmap = resList.getChi2CombinationCoeff2DMap(k, enablePlot=enablePlot, chi2dontloadThres=chi2dontloadThres,  zthres=zthres, opt_combination=opt_combination)
            coeffMap, sigma_coeff, b_coeff, xval, yval = resList.getBayesCombinationCoeff(indice=k, enablePlot=0, chi2dontloadThres=-1, zthres=zthres)
        except Exception as e:
            traceback.print_exc()
            if dont_skip_no_spc:
                stop
            else:
                continue

        idx_sigma = 2#coeffMap.shape[0]-1
        print("using idx_sigma = {}".format(idx_sigma))
        print("using sigma = {}".format(sigma_coeff[idx_sigma]))
        success = False
        success_end_reached = False
        for kr, c in enumerate(coeffMap[idx_sigma]):
            print("for kr={}, dtdnc={}, dtd={}, c={}".format(kr, xval, yval, c))
            point = [xval, b_coeff[kr]]
            #point = [yval, b_coeff[kr]]
            if c==1 and success==False:
                #point = [xval, b_coeff[k]*xval]
                print("adding point = {}".format(point))
                success_begin_point = point
                success = True
            elif c==0 and success==True:
                success_end_point = previous_point
                success_end_reached = True
                break
            
            previous_point = point
                
        if success==True and not success_end_reached:
            success_end_point = point
            
        if success == True:
            success_begin_points.append(success_begin_point)
            success_end_points.append(success_end_point)
            success_range.append([resList.list[k].name, xval, yval, success_begin_point[1], success_end_point[1]])
        else:
            success_range.append([resList.list[k].name, xval, yval, -1, -1])
            
        #check concordance between begin and end points
        if not len(success_begin_points) == len(success_end_points) :
            print("begin and end points have not the same dimension...")
            stooooop
            

        #print("coeff = {}".format(coeffMap))
        if coeffMapCumulative==None:
            coeffMapCumulative=np.zeros((coeffMap.shape[0],  coeffMap.shape[1]))
        
        coeffMapCumulative = np.add(coeffMapCumulative, coeffMap)
        plt.clf()
        plt.close()
        fig = plt.figure('dtreec {} coeff map'.format("sigma, b"), figsize=(9, 8))
        ax = fig.add_subplot(111)
        cmap = plt.get_cmap('RdYlGn') 
        ncolors = 200 #min(20,k+2)
        cmap = cmap_discretize(cmap, ncolors) 

        i = ax.matshow(np.transpose(coeffMapCumulative), interpolation='nearest', aspect='auto', cmap=cmap)
        
        xaxisticksinds = range(len(sigma_coeff))
        xaxisticks = ["{:.1f}".format(a) for a in sigma_coeff]
        plt.xticks(xaxisticksinds, xaxisticks, rotation='horizontal',verticalalignment='bottom') 
        
        yaxisticksinds = range(len(b_coeff))
        yaxisticks = ["{:.1e}".format(a) for a in b_coeff]
        plt.yticks(yaxisticksinds, yaxisticks, rotation='horizontal',verticalalignment='bottom')
        #plt.yticks(range(nlines), linesaxisticks, rotation='horizontal',verticalalignment='bottom')
        

        plt.xlabel('sigma coeff')
        plt.ylabel('b coeff')
        name1 = "Combination map - dataset optimized\nrelzErr-Threshold = {}\n(processed idx={}/{})".format(zthres, k+1, resList.n)
        plt.title(name1)
        #plt.legend(legendz)
        #plt.grid()
        #i = ax.imshow(image, interpolation='nearest')
        #if k==0:
        #    fig.colorbar(i)
        cbar = fig.colorbar(i)
        cbmax = np.amax(coeffMapCumulative)
        i.set_clim(vmin=cbmax-ncolors-0.5, vmax=cbmax+0.5)
        ncolorticksmaxcbar = 20+1
        cbarticks = np.linspace(cbmax-ncolors, cbmax, ncolorticksmaxcbar)
        cbar.set_ticks(cbarticks)
        cbar.set_ticklabels(["{:.1f}".format(a) for a in cbarticks])
        plt.draw()
        
        if enableExport:    
            tag = "dtreec_coeffoptim_b{}-sigma{:.2f}_relzerrthres{}".format("varying", sigma_coeff[idx_sigma], zthres)
            outdir = os.path.join(resList.analysisoutputdir, tag)
            if not os.path.exists(outdir):
                print("creating outputdir {}".format(outdir))
                os.makedirs(outdir)   
            outFigFile = os.path.join(outdir, 'coeffmapcumulative.png')
            plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
            
            outtxtFile = os.path.join(outdir, 'success_range_{}.txt'.format(tag))
            f = open(outtxtFile, 'w')
            for a in success_range:
                for b in a:
                    f.write("{}\t".format(b))
                f.write('\n')
            f.close()
            
            

        if enable_intermediate_plots:
            #NOTE: in order to deactivate showing this window at each iteration, use Agg (see imports)
            plt.clf()
            plt.close()
            
            fig = plt.figure('dtreec continuum sigma-coeff/b-coeff'.format(), figsize=(9, 8))
            ax = fig.add_subplot(111)
            ax.plot(b_coeff, np.transpose(coeffMap), label=b_coeff)     
            ax.xaxis.grid('True')
            ax.yaxis.grid('True')
            ax.yaxis.grid('True')
            ax.set_xscale('log')
            ax.set_ylim([-0.2, 1.2])
            #legend ?
            ax.set_title("spc name = {}".format(resList.list[k].name))
            plt.show()
            
        #display, plot coeff plots
        if 1:
            plt.clf()
            plt.close()
            
            fig = plt.figure('dtreec coeff points'.format(), figsize=(16, 8))
            ax = fig.add_subplot(111)
            xx = [a[0] for a in success_begin_points]
            yy = [a[1] for a in success_begin_points]
            ax.plot(xx, yy, 'xb')   
#            xx = [a[0] for a in failure_beginning_points]
#            yy = [a[1] for a in failure_beginning_points]
#            ax.plot(xx, yy, 'xr')   
            xx = [a[0] for a in success_end_points]
            yy = [a[1] for a in success_end_points]
            ax.plot(xx, yy, 'xr')  
            
            for a in range(len(success_begin_points)):
                x = success_begin_points[a][0]
                y1 = b_coeff[0]
                y2 = success_begin_points[a][1]
                z1 = success_end_points[a][1]
                z2 = b_coeff[len(b_coeff)-1]
                
                plt.plot((x, x), (y1,y2) , '-', color='lavender', label="{}_low".format(a) )
                plt.plot((x, x), (z1,z2) , '-', color='salmon', label="{}_up".format(a) )
            
            ax.xaxis.grid('True')
            ax.yaxis.grid('True')
            ax.yaxis.grid('True')
            #ax.set_xscale('log')
            #ax.set_ylim([-0.2, 1.2])
            ax.set_title("sigma_coeff={}\nspc in stats = {}/{}".format(sigma_coeff[idx_sigma], len(success_begin_points), k+1))
            plt.draw()
            if enableExport:    
                outdir = os.path.join(resList.analysisoutputdir, tag)
                if not os.path.exists(outdir):
                    print("creating outputdir {}".format(outdir))
                    os.makedirs(outdir)   
                outFigFile = os.path.join(outdir, 'coeffpoints.png')
                plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
                
            
            outtxtFile = os.path.join(outdir, 'success_range_{}.txt'.format(tag))
            f = open(outtxtFile, 'w')
            for a in success_range:
                for b in a:
                    f.write("{}\t".format(b))
                f.write('\n')
            f.close()
            
        time.sleep(0.05)
        
    #plot in 3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    _colors = ['b', 'r']
    _markers = ['^', 'o']
    x_begin = []
    y_begin = []
    z_begin = []
    x_end = []
    y_end = []
    z_end = []
    for a in range(len(success_range)):
        x = success_range[a][1]
        y = success_range[a][2]
        z = success_range[a][3]
        x_begin.append(x)
        y_begin.append(y)
        z_begin.append(z)
        z = success_range[a][4]
        x_end.append(x)
        y_end.append(y)
        z_end.append(z)
    ax.scatter(x_begin, y_begin, z_begin, c=_colors[0], marker=_markers[0])
    ax.scatter(x_end, y_end, z_end, c=_colors[1], marker=_markers[1])
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    plt.show()
        
        
def exportAmazed03SVMTable(resDir, diffthres, spcName=""):
    print('using amazed results full path: {0}'.format(resDir))
    zrefmin=0.01
    zrefmax=10.0
    
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName, methodName="", zrefmin=zrefmin, zrefmax=zrefmax)
    if resList.n <1:
        print('No results loaded...')
        return

    respa = rp.ResParser(resDir) 
    idxFlag = 0
    idxSpcName = 1
    idxRedshiftCandidate = 2
    idxRelRedshiftErrBestCandidate = 3
    idxMeritCandidate_lm = 4
    idxRelMeritBestCandidate_lm = 5
    idxMeritStd = 6
    idxNCandSignificant = 7
    
    
    idxMeritCandidate_chi2= 8
    idxRelMeritBestCandidate_chi2= 9
    
    idxDTDNoContinuum = 10
    idxDTD = 11
    idxDTDNC_over_idxDTD = 12
    
    N = 13
    data_list = []
    
    #nres = 100
    nres = resList.n
    nSkippedNoZref = 0
    nSkippedNoDTD = 0
    for k in range(nres):
        print("\nprocessing result #{}/{}".format(k+1, resList.n))
        zref = resList.list[k].zref
        print("zref is {}".format(zref))
        
        #### get the LM and Chi2continuum candidates
        candid_zvalCandidates, candid_displayParamsBundle = respa.getAutoCandidatesList(resList.list[k].name)
        lm_candid_displayParamsBundleIdx = 0
        chi2_candid_displayParamsBundleIdx = 1

        #get dtranspose
        filepaths, filenames = respa.getAutoChi2FullPath(resList.list[k].name)
        filepath = filepaths[1] #1=linemodel
        print("Found chi2 filepath: {}".format(filepath))
        if not os.path.exists(filepath):
            print("Problem while retrieving chi2 filepath.. using: {}".format(filepath))
            stoooooop
        else:
            print("using Chi2 file path : ".format(filepath))
            chi2FileResult = chisq.ResultChisquare(filepath)
            dtransposeD = chi2FileResult.amazed_dTransposeD[0]
            dtransposeDNoContinuum = chi2FileResult.amazed_dTransposeDNoContinuum[0] #chi2.getFluxMedian()
            if dtransposeD==0 or dtransposeDNoContinuum==0:
                print("Problem while retrieving dtransposeD: {}".format(dtransposeD))
                nSkippedNoDTD += 1
                continue
                
            
        #converting to numpy arrays
        z_lm = np.array(candid_zvalCandidates)
        nz = len(z_lm)
        #print("z_lm = {}".format(z_lm))
        chi2_lm = np.array(candid_displayParamsBundle[lm_candid_displayParamsBundleIdx]['merits'])
        #print("chi2_lm = {}".format(chi2_lm))
        chi2_continuum = np.array(candid_displayParamsBundle[chi2_candid_displayParamsBundleIdx]['merits'])
        #print("chi2_continuum = {}".format(chi2_continuum))

        #print        
        for i, a in enumerate(z_lm):
            print("#{} : redshifts (z_lm)={:<15}, chi2_lm={:<15}, chi2_continuum={:<15}".format(i, z_lm[i], chi2_lm[i], chi2_continuum[i]))       
        #######
              
              
        #get the ref candidate              
        print("Searching for candidate with zref = {}".format(zref))
        thres_zref_zerr = 1e-2
        indsZrefExtremum = [i for i,z in enumerate(z_lm) if abs(z-zref)<thres_zref_zerr]
        if len(indsZrefExtremum)!=1: #if there is more than 1 correct candidate for that source, skip, somethind is wrong...
            print("ERROR: found more than 1 zref candidate...")
            nSkippedNoZref += 1
            continue
        iZrefExtremum = indsZrefExtremum[0]
        print("Found zref candidate idx = {}".format(iZrefExtremum))

        #get the wrong candidates        
        thres_zwrong_zerr_min = 0.001*(1+zref)
        thres_zwrong_zerr_max = 10.0*(1+zref)
        thres_extrema_id = 10
        indsZwrongExtrema = [i for i,z in enumerate(z_lm) if abs(z-zref)>thres_zwrong_zerr_min and abs(z-zref)<thres_zwrong_zerr_max and i<thres_extrema_id]
   
        print("inds extrema for zref found = {}".format(indsZrefExtremum))
        print("inds extrema for zwrong found = {}".format(indsZwrongExtrema))
        
        #modify merits
        meritsModified_lm = [(m-dtransposeDNoContinuum)/dtransposeDNoContinuum for m in chi2_lm]
        meritsModified_chi2 = [(m-dtransposeD)/dtransposeD for m in chi2_continuum]         
         
        #get the significant condidates
        thres_significant = 0.2
        indsSignificant = [i for i, m in enumerate(meritsModified_lm) if m/meritsModified_lm[iZrefExtremum]>thres_significant]
        indBestMerit = 0  
        
        #build the common (all candidates for this source) template vector
        _data = np.zeros((N))
        _data[idxSpcName] = resList.list[k].name.split("_")[1]
        _data[idxMeritStd] = chi2FileResult.getFluxStd()/dtransposeDNoContinuum
        _data[idxNCandSignificant] = len(indsSignificant)
        _data[idxDTDNoContinuum] = dtransposeDNoContinuum
        _data[idxDTD] = dtransposeD
        _data[idxDTDNC_over_idxDTD] = dtransposeDNoContinuum/dtransposeD
        

        #build the data
        for i, kInd in enumerate(indsSignificant):
            data = np.copy(_data) 
            if kInd == iZrefExtremum:
                data[idxFlag] = 1
            else:
                data[idxFlag] = 0
                
            data[idxRedshiftCandidate] = z_lm[kInd]
            data[idxRelRedshiftErrBestCandidate] = (z_lm[kInd]-z_lm[indBestMerit])/(1+z_lm[indBestMerit])
            
            data[idxMeritCandidate_lm] = meritsModified_lm[kInd]
            data[idxRelMeritBestCandidate_lm] = meritsModified_lm[kInd]-meritsModified_lm[indBestMerit]
            
            data[idxMeritCandidate_chi2] = meritsModified_chi2[kInd]
            data[idxRelMeritBestCandidate_chi2] = meritsModified_chi2[kInd]-meritsModified_chi2[indBestMerit]


                        
            data_list.append(data)
        
        print("data = {}".format(data_list))
    
    print("Finished extraction (nres={}) : nSkippedNoZref={}, nSkippedNoDTD={}".format(nres, nSkippedNoZref, nSkippedNoDTD))
    enableExport = 1 
    if enableExport:
        outdir = os.path.join(resList.analysisoutputdir, "amazed03_svm")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        outFile = os.path.join(outdir, 'svm_data.txt'.format())
        np.savetxt(outFile, data_list)  
         
       
   
def plotReducedZcandidates(resDir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=0.01, opt='brief')
    
    resList.plotReducedZcandidates("amazed", True)     
        
def exportAnalysis(resDir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1, opt='full')
    nextrema = 5;
    #0 names
    names = [a.name for a in resList.list]
    #1 zref
    zref = [a.zref for a in resList.list]
    #2 zcalc
    zcalc = [a.zcalc for a in resList.list]
    #3 merit
    meritCalc = [a.chi2 for a in resList.list]
    #4 tpl
    tplZCalc = resList.getBestTplZCalc()
    #5
    meritRef = resList.getBestChi2ZRef()
    #6
    tplZRef = resList.getBestTplZRef()
    #7 
    zrefFoundChi2Raw = resList.getClosestZcandidatesZrefList("raw", "amazed", enablePlot=False, enableExport=True, nextrema = nextrema)
    #8
    zrefFoundChi2NoContinuum = resList.getClosestZcandidatesZrefList("nocontinuum", "amazed", enablePlot=False, enableExport=True, nextrema = nextrema)
    #9
    zrefFoundCorr = resList.getClosestZcandidatesZrefList("corr", "amazed", enablePlot=False, enableExport=True, nextrema = nextrema)
    #10
    tplMissingRate = resList.tplMissingRate_to_liste()
    #11
    relposraw = resList.getClosestZcandidateRelativePositionList("raw", "amazed", enablePlot=False, enableExport=True, nextrema = nextrema)
    #12
    relposnocont = resList.getClosestZcandidateRelativePositionList("nocontinuum", "amazed", enablePlot=False, enableExport=True, nextrema = nextrema)
    #13
    relposcorr = resList.getClosestZcandidateRelativePositionList("corr", "amazed", enablePlot=False, enableExport=True, nextrema = nextrema)
    
    outFile = os.path.join(resList.analysisoutputdir, "analysis_nextrema{}.csv".format(nextrema))
    fout = open( outFile, "w" ) 
    fout.write( "name\tzref\tzcalc\
    \tmeritZCalc\ttpl for zcalc\
    \tmeritZRef\ttpl for zref\
    \tzclosest in chi2raw\tzclosest in chi2noncontinuum\tzclosest in corr\
    \ttpl missingrate\
    \tzrelpos in chi2raw\tzrelpos in chi2noncontinuum\tzrelpos in corr\n" )
    for k in range(resList.n):
        fout.write( "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n".format(\
        names[k], zref[k], zcalc[k]\
        , meritCalc[k], tplZCalc[k]\
        , meritRef[k],tplZRef[k]\
        , zrefFoundChi2Raw[k], zrefFoundChi2NoContinuum[k], zrefFoundCorr[k]\
        , tplMissingRate[k]\
        , relposraw[k],  relposnocont[k],  relposcorr[k]))
    fout.close()

def plotClosestZCandidates(resDir, diffthres, nextrema=5, extremaType = "linemodel"):    
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief')
    zrefFoundChi2Raw = resList.getClosestZcandidatesZrefList(extremaType, "amazed", enablePlot=True, enableExport=True, nextrema = nextrema)
    #zrefFoundChi2Raw = resList.getClosestZcandidatesZrefList("raw", "", True)

def plotPrecisionHist(resDir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1.0, opt='brief')
    
    zrefMin = 0.0
    zrefMax = 20.0
    relzerrList = []
    
    nres = resList.n
    for k in range(nres):
        print("\nprocessing result #{}/{}".format(k+1, resList.n))
        zref = resList.list[k].zref
        zcalc = resList.list[k].zcalc
        relzerr = (resList.list[k].zref-zcalc)/(1.0+resList.list[k].zref)
        if zref >= zrefMin and zref <= zrefMax:
            if abs(relzerr) < 5e-3:
                relzerrList.append(relzerr)
    
    plt.figure("Precision Histogram", figsize=(10,8))
    #xvect = range(len(relzerrList))
    nbins = 40
    plt.hist(relzerrList, bins=nbins, normed=0, histtype='step')
    plt.xlim([-5e-3, 5e-3])
    plt.grid(True) # Affiche la grille
    plt.xlabel('(zcalc-zref)/(1+zref)')
    plt.ylabel('count')
    plt.title('Histogram\n{}/{} included results'.format(len(relzerrList), nres))
    plt.show()
    
def plotVelocityError(resDir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1.0, opt='brief')
    
    elvList = []
    alvList = []
    idXExtremum = 0
    elvErrorList = []
    elvRefList = []
    
    zrefMin = 0.0
    zrefMax = 20.0
    relzerrList = []
    relzerr_threshold = 100#5e-3
    
    nres = resList.n
    for k in range(nres):
        print("\nprocessing result #{}/{}".format(k+1, resList.n))
        zref = resList.list[k].zref
        zcalc = resList.list[k].zcalc
        relzerr = (resList.list[k].zref-zcalc)/(1.0+resList.list[k].zref)
        if zref >= zrefMin and zref <= zrefMax:
            if abs(relzerr) < relzerr_threshold:
                relzerrList.append(abs(relzerr))
                lmModelFitPath = resList.resParser.getLineModelResultPath(resList.list[k].name, idXExtremum)
                if os.path.exists(lmModelFitPath):
                    mres = modelresult.ModelResult(lmModelFitPath)
                else:
                    print("ERROR: unable to find the modelResult file... aborting...")
                    stop
                try:
                    elv = mres.vel_emission
                    alv = mres.vel_absorption
                except:
                    print("ERROR-EXCEPT: unable to find the line velocities file... aborting...")
                    stop
                elvList.append(elv) 
                alvList.append(alv)
                elvRefList.append(resList.list[k].refValues['elvelocity'])
                elvErrorList.append((elv-resList.list[k].refValues['elvelocity']))

    print("got relzerr n = {}".format(len(relzerrList)))
    print("got elvError n = {}".format(len(elvErrorList)))        
      
    enableExport = 1 
    enablePlot = 1
    if enableExport or enablePlot:
        outdir = os.path.join(resList.analysisoutputdir, "velocity")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
            plt.figure("Velocity Error Histogram", figsize=(10,8))
            
        #plotting Errorvelocity histogram
        nbins = 80
        plt.hist(elvErrorList, bins=nbins, normed=0, histtype='stepfilled')
        #plt.xlim([-5e-3, 5e-3])
        plt.grid(True) # Affiche la grille
        plt.xlabel('abs(ELVcalc-ELVref)')
        plt.ylabel('count')
        plt.title('Histogram\n{}/{} included results'.format(len(elvErrorList), nres))
        if enableExport:
            outFigFile = os.path.join(outdir, 'velocityErr_hist.png'.format())
            plt.savefig( outFigFile, bbox_inches='tight')
        if enablePlot:
            plt.show() 
            
        #Plotting hist versus velocity
        #elvAbsErrorList = [np.abs(a) for a in elvErrorList]
        outFileNoExt = 'stats_versusSigma_hist' 
        outFilepathNoExt = os.path.join(outdir,outFileNoExt)
        lstats.PlotAmazedVersusBinsHistogram(relzerrList, elvErrorList, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=enableExport, mtype='ERROR_SIGMA', nPercentileDepth=2) 
    
  
def plotTplMissingRate(resDir, opt=0):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=0.01)
    #resList.chi2_export_liste()
    resList.tplMissingRate_plot()      
        
def plotChi2XYPlane(resDir, opt=0):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir)
    #resList.chi2_export_liste()
    resList.chi2_plot()  
  
def comparePeakDetection(resDir, refresdir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1, opt='brief')
    refresList = ResultList(refresdir, diffthreshold=-1, opt='brief')
    resList.getPeakDetectionComparison(refresList) 
 
 
def exportLineDetectionStats(resDir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1, opt='brief')
    resList.getLineDetectionStatsList()
      
def compareLineDetection(resDir, refresdir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1, opt='brief')
    refresList = ResultList(refresdir, diffthreshold=-1, opt='brief')
    resList.getLineDetectionComparison(refresList)
  
    
def exportLineMatchingStats(resDir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=-1, opt='brief')
    resList.getLineMatchingStatsList(enablePlot = False)
   
def compareFailures(resDir, refresdir, diffthreshold=0.01, zrefmin=-1, zrefmax=20):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthreshold, opt='brief', spcName="", methodName="", zrefmin=zrefmin, zrefmax=zrefmax)
    refresList = ResultList(refresdir, diffthreshold=diffthreshold, opt='brief', spcName="", methodName="", zrefmin=zrefmin, zrefmax=zrefmax)
    resList.getFailuresComparison(refresList)
  
#def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
#    xo = float(xo)
#    yo = float(yo)    
#    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
#    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
#    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
#    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
#                            + c*((y-yo)**2)))
#    return g.ravel()
    
def plotContinuumIndexes(resDir, diffthres, spcName=""):
    print('using amazed results full path: {0}'.format(resDir))
    print('using spc filter by name: {}'.format(spcName))
    zrefmin=0.01
    zrefmax=1.5
    zrefmin=-1
    zrefmax=50.0
    
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName, methodName="", zrefmin=zrefmin, zrefmax=zrefmax)
    if resList.n <1:
        print('No results loaded...')
        return

    s = rp.ResParser(resDir)        
    continuumIndexesZrefList = []   
    zrefZrefList = []  
    continuumIndexesZwrongList = []
    zrefZwrongList = []
    
    #if processAmazedatZ has been used, both of these has to be 1
    thres_extrema_id = 1#10
    thres_extremasearch_id = 1#10
    
    nSkippedZref = 0
    #nres = 20
    nres = resList.n
    for k in range(nres):
        print("\nprocessing result #{}/{}".format(k+1, resList.n))
        zref = resList.list[k].zref
        print("zref is {}".format(zref))
        redshifts, merits = resList.getZCandidatesFromAmazedChi2Extrema(k, chi2Type="linemodel", nextrema=thres_extremasearch_id, enableZrangePerTpl=False)
        print("redshifts: {}".format(redshifts))       
        print("merits: {}".format(merits))  
        thres_zref_zerr = 1e-2
        indsZrefExtremum = [i for i,z in enumerate(redshifts) if abs(z-zref)<thres_zref_zerr]
        if 0: #keep only the potential Ha/OII ambiguities
            thres_zwrong_zerr_min = 0.05*(1+zref)
            thres_zwrong_zerr_max = 1.5*(1+zref)
        else:
            thres_zwrong_zerr_min = 0.001*(1+zref)
            thres_zwrong_zerr_max = 10.0*(1+zref)        
        indsZwrongExtrema = [i for i,z in enumerate(redshifts) if abs(z-zref)>thres_zwrong_zerr_min and abs(z-zref)<thres_zwrong_zerr_max and i<thres_extrema_id]
        
        print("inds extrema for zref found = {}".format(indsZrefExtremum))
        print("inds extrema for zwrong found = {}".format(indsZwrongExtrema))
        if len(indsZrefExtremum)>0 or len(indsZwrongExtrema)>0:
            filepaths, filenames = s.getAutoChi2FullPath(resList.list[k].name)
            #filepath = filepaths[0] #0 for linemodel method
            filepath = filepaths[1] #1 for amazed03 method
            
            print("Found chi2 filepath: {}".format(filepath))
            if not os.path.exists(filepath):
                print("Problem while retrieving chi2 filepath.. using: {}".format(filepath))
                continue
            else:
                print("using Chi2 file path : ".format(filepath))
            chi2 = chisq.ResultChisquare(filepath)
            #print("chi2 is: {}".format(chi2))
            
        if len(indsZrefExtremum)>0:
            for ie in indsZrefExtremum:
                continuumIndexesZref = chi2.amazed_continuumIndexes[ie]
                print("continuum indexes zref are: {}".format(continuumIndexesZref))
                continuumIndexesZrefList.append(continuumIndexesZref)
                zrefZrefList.append(zref)
        else:
            nSkippedZref+=1
            
        if len(indsZwrongExtrema)>0:
            for ie in indsZwrongExtrema:
                continuumIndexesZwrong = chi2.amazed_continuumIndexes[ie]
                print("continuum indexes zwrong are: {}".format(continuumIndexesZwrong))
                continuumIndexesZwrongList.append(continuumIndexesZwrong)
                zrefZwrongList.append(zref)
            
    print("\nINFO-SKIPPED: nskipped zref = {}".format(nSkippedZref))
    
    
    enablePlot = 0
    enableExport = 1
    enablePlotZref = 1
    enablePlotZwrong = 1
    enablePlotLineSeparation = 1
    
    if enablePlot or enableExport:
        
        idx_continuum_index_Lya = 0
        idx_continuum_index_OII = 1
        idx_continuum_index_OIII = 2
        idx_continuum_index_Ha = 3
        idx_continuum_index_CIV = 4
        idx_continuum_index_CIII = 5
        
        fig = plt.figure('continuum indexes - zref'.format(), figsize=(16,12))
        ax = fig.add_subplot(111)
        
        figWrong = plt.figure('continuum indexes - zwrong'.format(), figsize=(16,12))
        axWrong = figWrong.add_subplot(111)
        
        if enablePlotZref:
            #xvect = range(len(continuumIndexesZrefList))
            #xvect = zrefZrefList
            #print continuumIndexesZrefList[0]['color']
            xvect = [a['break'][idx_continuum_index_Lya] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_Lya]) and not np.isnan(a['color'][idx_continuum_index_Lya])]
            yvect = [a['color'][idx_continuum_index_Lya] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_Lya]) and not np.isnan(a['color'][idx_continuum_index_Lya])]
            print("1. xvect={}".format(xvect))
            print("1. yvect={}".format(yvect))
            ax.plot(xvect, yvect, 'xk', label='zref Lya')
            
            xvect = [a['break'][idx_continuum_index_OII] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_OII]) and not np.isnan(a['color'][idx_continuum_index_OII])]
            yvect = [a['color'][idx_continuum_index_OII] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_OII]) and not np.isnan(a['color'][idx_continuum_index_OII])]
            print("1. xvect={}".format(xvect))
            print("1. yvect={}".format(yvect))
            ax.plot(xvect, yvect, 'xb', label='zref OII')
            
            xvect = [a['break'][idx_continuum_index_OIII] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_OIII]) and not np.isnan(a['color'][idx_continuum_index_OIII])]
            yvect = [a['color'][idx_continuum_index_OIII] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_OIII]) and not np.isnan(a['color'][idx_continuum_index_OIII])]
            print("2. xvect={}".format(xvect))
            print("2. yvect={}".format(yvect))
            ax.plot(xvect, yvect, 'xg', label='zref OIII')
            
            xvect = [a['break'][idx_continuum_index_Ha] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_Ha]) and not np.isnan(a['color'][idx_continuum_index_Ha])]
            yvect = [a['color'][idx_continuum_index_Ha] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_Ha]) and not np.isnan(a['color'][idx_continuum_index_Ha])]
            print("3. xvect={}".format(xvect))
            print("3. yvect={}".format(yvect))
            ax.plot(xvect, yvect, 'xr', label='zref Ha')
            
            xvect = [a['break'][idx_continuum_index_CIV] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_CIV]) and not np.isnan(a['color'][idx_continuum_index_CIV])]
            yvect = [a['color'][idx_continuum_index_CIV] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_CIV]) and not np.isnan(a['color'][idx_continuum_index_CIV])]
            print("4. xvect={}".format(xvect))
            print("4. yvect={}".format(yvect))
            ax.plot(xvect, yvect, 'xm', label='zref CIV')
            
            xvect = [a['break'][idx_continuum_index_CIII] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_CIII]) and not np.isnan(a['color'][idx_continuum_index_CIII])]
            yvect = [a['color'][idx_continuum_index_CIII] for a in continuumIndexesZrefList if not np.isnan(a['break'][idx_continuum_index_CIII]) and not np.isnan(a['color'][idx_continuum_index_CIII])]
            print("5. xvect={}".format(xvect))
            print("5. yvect={}".format(yvect))
            ax.plot(xvect, yvect, 'xy', label='zref CIII')
            
            ax.xaxis.grid(True,'major')
            ax.yaxis.grid(True,'major')
            ax.legend()
            ax.set_ylabel('color')
            ax.set_xlabel('break')
            titleStr = "Continuum Indexes - zref, nSkippedZref = {}/{}".format(nSkippedZref, nres)
            ax.set_title(titleStr)
            
            
        if enablePlotZwrong:
            #xvect = range(len(continuumIndexesZwrongList))
            #xvect = zrefZwrongList
            xvect = [a['break'][idx_continuum_index_OII] for a in continuumIndexesZwrongList if not np.isnan(a['break'][idx_continuum_index_OII]) and not np.isnan(a['color'][idx_continuum_index_OII])]
            yvect = [a['color'][idx_continuum_index_OII] for a in continuumIndexesZwrongList if not np.isnan(a['break'][idx_continuum_index_OII]) and not np.isnan(a['color'][idx_continuum_index_OII])]
            print("4. xvect={}".format(xvect))
            print("4. yvect={}".format(yvect))
            axWrong.plot(xvect, yvect, 'vb', label='zwrong OII')
            
            
            xvect = [a['break'][idx_continuum_index_OIII] for a in continuumIndexesZwrongList if not np.isnan(a['break'][idx_continuum_index_OIII]) and not np.isnan(a['color'][idx_continuum_index_OIII])]
            yvect = [a['color'][idx_continuum_index_OIII] for a in continuumIndexesZwrongList if not np.isnan(a['break'][idx_continuum_index_OIII]) and not np.isnan(a['color'][idx_continuum_index_OIII])]
            print("5. xvect={}".format(xvect))
            print("5. yvect={}".format(yvect))
            axWrong.plot(xvect, yvect, 'vg', label='zwrong OIII')
            
            xvect = [a['break'][idx_continuum_index_Ha] for a in continuumIndexesZwrongList if not np.isnan(a['break'][idx_continuum_index_Ha]) and not np.isnan(a['color'][idx_continuum_index_Ha])]
            yvect = [a['color'][idx_continuum_index_Ha] for a in continuumIndexesZwrongList if not np.isnan(a['break'][idx_continuum_index_Ha]) and not np.isnan(a['color'][idx_continuum_index_Ha])]
            print("6. xvect={}".format(xvect))
            print("6. yvect={}".format(yvect))
            axWrong.plot(xvect, yvect, 'vr',label='zwrong Ha')
            
            axWrong.xaxis.grid(True,'major')
            axWrong.yaxis.grid(True,'major')
            axWrong.legend()
            axWrong.set_ylabel('color')
            axWrong.set_xlabel('break')
            titleStr = "Continuum Indexes - zwrong, nSkippedZref = {}/{}".format(nSkippedZref, nres)
            axWrong.set_title(titleStr)
            
        if enablePlotLineSeparation:
            a = 1.0
            b = 0.25
            xvect = [-1.0, 0.0, 1.0]
            yvect = [a*x+b for x in xvect]
            if enablePlotZref:
                ax.plot(xvect, yvect, 'k-', linestyle = "dashed", label='line: a={}, b={}'.format(a, b))
            if enablePlotZref:
                axWrong.plot(xvect, yvect, 'k-', linestyle = "dashed", label='line: a={}, b={}'.format(a, b))
        
        #plt.xlim([1e-5, 10])
        #plt.ylim([0, 100])
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        #ax.set_xscale('log')
        #plt.grid(True) # Affiche la grille
        #plt.legend()
        #plt.ylabel('color')
        #plt.xlabel('break')
        #titleStr = "Continuum Indexes, nSkippedZref = {}/{}".format(nSkippedZref, nres)
        #plt.title(titleStr)
        
        if enableExport:
            outdir = os.path.join(resList.analysisoutputdir, "continuum_indexes")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
                
            if enablePlotZref:
                outFigFile = os.path.join(outdir, 'continuum_indexes_zref.png'.format())
                fig.savefig( outFigFile, bbox_inches='tight')
            if enablePlotZwrong:
                outFigFile = os.path.join(outdir, 'continuum_indexes_zwrong.png'.format())
                figWrong.savefig( outFigFile, bbox_inches='tight')
            
            #export all continuum indexes in files
            for k in range(6):
                fPath = os.path.join(outdir, "scatter_ycolor_xbreak_{}.csv".format(k))
                xvect = [a['break'][k] for a in continuumIndexesZrefList if not np.isnan(a['break'][k]) and not np.isnan(a['color'][k])]
                yvect = [a['color'][k] for a in continuumIndexesZrefList if not np.isnan(a['break'][k]) and not np.isnan(a['color'][k])]
                f = open(fPath, 'w')
                for kx, x in enumerate(xvect):
                    f.write("{}\t{}\n".format(xvect[kx], yvect[kx]))
                f.close()
                
                #generate heat map
                xvect.append(-4)
                yvect.append(-2)
                xvect.append(-4)
                yvect.append(3)
                xvect.append(2)
                yvect.append(-2)
                heatmap, xedges, yedges = np.histogram2d(xvect, yvect, bins=50)
                extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                
                heatmap_blurred = ndimage.gaussian_filter(heatmap, sigma=2)
                maxi = np.max(heatmap_blurred)
                print("maxi = {}".format(maxi))
                heatmap_blurred = heatmap_blurred/maxi
                #interp on reg fixed grid
                step = 0.1
                xmin = -3
                xmax = 2
                n = (xmax-xmin)/step
                grid_x, grid_y = np.mgrid[xmin: xmax: step, xmin: xmax: step]
                print("shape: gridx={}, gridy={}".format(grid_x.shape, grid_y.shape))
                points = []
                values = []
                n1 = len(heatmap_blurred)
                n2 = len(heatmap_blurred[0])
                step1 = xedges[1]-xedges[0]
                step2 = yedges[1]-yedges[0]
                for k1 in range(n1):
                    for k2 in range(n2):
                        points.append((k1*step1+xedges[0], k2*step2+yedges[0]))
                        values.append(heatmap_blurred[k1][k2])
                heatmap_blurred_fixedtbl = griddata(points, values, (grid_x, grid_y), method='linear', fill_value=0.0)
                extentfixed = [xmin, xmax, xmin, xmax]
                        
                #plt.figure()
                plt.matshow(heatmap.T, extent=extent, interpolation='nearest', origin='lower')
                plt.ylabel('color')
                plt.xlabel('break')
                plt.colorbar()
                outFigFileHmap = os.path.join(outdir, "heatmap_ycolor_xbreak_{}.png".format(k))
                plt.savefig( outFigFileHmap, bbox_inches='tight')
                
                plt.matshow(heatmap_blurred.T, extent=extent, interpolation='nearest', origin='lower')
                plt.ylabel('color')
                plt.xlabel('break')
                plt.colorbar()
                outFigFileHmap = os.path.join(outdir, "heatmap_ycolor_xbreak__blurred_{}.png".format(k))
                plt.savefig( outFigFileHmap, bbox_inches='tight')
                #plt.show()
                
                plt.matshow(heatmap_blurred_fixedtbl.T, extent=extentfixed, interpolation='nearest', origin='lower')
                plt.ylabel('color')
                plt.xlabel('break')
                plt.colorbar()
                outFigFileHmap = os.path.join(outdir, "heatmap_ycolor_xbreak__blurred_fixedtbl_{}.png".format(k))
                plt.savefig( outFigFileHmap, bbox_inches='tight')
                #plt.show()

                outFileConfigHmap = os.path.join(outdir, "heatmap_ycolor_xbreak__blurred_fixedgrid_{}.cfg".format(k)) 
                f = open(outFileConfigHmap, 'w')
                f.write("{}\t{}\t{}".format(xmin, xmax, step))
                f.write("\n")
                f.write("{}\t{}\t{}".format(xmin, xmax, step))
                f.close()
                
                outFiledataHmap = os.path.join(outdir, "heatmap_ycolor_xbreak__blurred_fixedgrid_{}.dat".format(k)) 
                f = open(outFiledataHmap, 'w')
                for xd in range(len(heatmap_blurred_fixedtbl)):
                    for yd in range(len(heatmap_blurred_fixedtbl[0])):
                        f.write("{}\t".format(heatmap_blurred_fixedtbl[xd][yd]))
                    f.write("\n")
                f.close()
                
#                #fit
#                # add some noise to the data and try to fit the data generated beforehand
#                #initial_guess = (10,0,0,1,1,0,1)
#                #popt, pcov = opt.curve_fit(twoD_Gaussian, (xedges, yedges), data_noisy, p0=initial_guess)
#                if 0:                
#                    #create data
#                    x = np.linspace(0, 200, 201)
#                    y = np.linspace(0, 200, 201)
#                    x, y = np.meshgrid(x, y)
#                    data = twoD_Gaussian((x, y), 10.0,0.1,0.1,1.0,1.0,0.1,1.0)
#                    # plot twoD_Gaussian data generated above
#                    #plt.figure()
#                    plt.imshow(data.reshape(201, 201))
#                    plt.colorbar()

            
                
            
        if enablePlot:
            plt.show() 

def exportNLinesExternal(resDir, diffthres, spcName=""):
    print('using amazed results full path: {0}'.format(resDir))
    print('using spc filter by name: {}'.format(spcName))
    
    
    zrefmin=-1
    zrefmax=50.0
    
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName, methodName="", zrefmin=zrefmin, zrefmax=zrefmax)
    if resList.n <1:
        print('No results loaded...')
        return
       
    nlinesZrefList = []  
    
    thres_extremasearch_id = 10
    
    nSkippedZref = 0
    #nres = 20
    nres = resList.n
    for k in range(nres):
        print("\nprocessing result #{}/{}".format(k+1, resList.n))
        zref = resList.list[k].zref
        print("zref is {}".format(zref))
        redshifts, merits = resList.getZCandidatesFromAmazedChi2Extrema(k, chi2Type="linemodel", nextrema=thres_extremasearch_id, enableZrangePerTpl=False)
        print("redshifts: {}".format(redshifts))       
        print("merits: {}".format(merits))  
        thres_zref_zerr = 1e-2
        indsZrefExtremum = [i for i,z in enumerate(redshifts) if abs(z-zref)<thres_zref_zerr]
        
        print("inds extrema for zref found = {}".format(indsZrefExtremum))
            
        if len(indsZrefExtremum)>0:
            ie = indsZrefExtremum[0]
            lmModelFitPath = resList.resParser.getLineModelResultPath(resList.list[k].name, ie)
            if os.path.exists(lmModelFitPath):
                mres = modelresult.ModelResult(lmModelFitPath)
            else:
                print("ERROR: unable to find the modelResult file... aborting...")
                nSkippedZref+=1
                nlinesZrefList.append(-1)
                continue
            res_nlines = mres.getNLinesStrong(redshifts[ie])
                
            nlinesZrefList.append(res_nlines)
        else:
            nSkippedZref+=1
            nlinesZrefList.append(-1)
            continue
            
    print("\nINFO-SKIPPED: nskipped zref = {}".format(nSkippedZref))
    
    enableExport = 1 
    if enableExport:
        outdir = os.path.join(resList.analysisoutputdir, "lm_nlinesvalid")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
                        
        outFileNoExt = 'external_lmnlinesstrong' 
        outFilepathNoExt = os.path.join(outdir,outFileNoExt)
        outFile = '{}.txt'.format(outFilepathNoExt)
        f = open(outFile, "w")
        for k in range(nres):
            f.write("{}\t{}\n".format(resList.list[k].name, nlinesZrefList[k]))
        f.close()    
    
def exportSVMTable(resDir, diffthres, spcName=""):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName, methodName="", zrefmin=0.01, zrefmax=1.5)
    if resList.n <1:
        print('No results loaded...')
        return

    resultParser = rp.ResParser(resDir) 
    idxFlag = 0
    idxSpcName = 1
    idxRedshiftCandidate = 2
    idxRelRedshiftErrBestCandidate = 3
    idxMeritCandidate = 4
    idxRelMeritBestCandidate = 5
    idxMeritStd = 6
    idxNCandSignificant = 7
    
    #idxContIndexColorLya = 12
    #idxContIndexBreakLya = 13
    
    idxContIndexColorOII = 8
    idxContIndexBreakOII = 9
    idxContIndexColorOIII = 10
    idxContIndexBreakOIII = 11
    idxContIndexColorHa = 12
    idxContIndexBreakHa = 13
    
    N = 14
    data_list = []
    
    #nres = 100
    nres = resList.n
    for k in range(nres):
        print("\nprocessing result #{}/{}".format(k+1, resList.n))
        zref = resList.list[k].zref
        print("zref is {}".format(zref))
        redshifts, merits = resList.getZCandidatesFromAmazedChi2Extrema(k, chi2Type="linemodel", nextrema=10, enableZrangePerTpl=False)
        print("redshifts: {}".format(redshifts))       
        print("merits: {}".format(merits))  
        
        #get the ref candidate
        thres_zref_zerr = 1e-1
        indsZrefExtremum = [i for i,z in enumerate(redshifts) if abs(z-zref)<thres_zref_zerr]
        if len(indsZrefExtremum)!=1: #if there is more than 1 correct candidate for that source, skip, somethind is wrong...
            continue
        iZrefExtremum = indsZrefExtremum[0]

        #get the wrong candidates        
        thres_zwrong_zerr_min = 0.001*(1+zref)
        thres_zwrong_zerr_max = 10.0*(1+zref)
        thres_extrema_id = 5
        indsZwrongExtrema = [i for i,z in enumerate(redshifts) if abs(z-zref)>thres_zwrong_zerr_min and abs(z-zref)<thres_zwrong_zerr_max and i<thres_extrema_id]
   
        print("inds extrema for zref found = {}".format(indsZrefExtremum))
        print("inds extrema for zwrong found = {}".format(indsZwrongExtrema))
        
        filepaths, filenames = resultParser.getAutoChi2FullPath(resList.list[k].name)
        filepath = filepaths[0]
        print("Found chi2 filepath: {}".format(filepath))
        if not os.path.exists(filepath):
            print("Problem while retrieving chi2 filepath.. using: {}".format(filepath))
            continue
        else:
            print("using Chi2 file path : ".format(filepath))
        chi2 = chisq.ResultChisquare(filepath)
        #print("chi2 is: {}".format(chi2))
        chi2Std = chi2.getFluxStd()
        chi2Med = chi2.getFluxMedian()
        #modify merits
        meritsModified = [m-chi2Med for m in merits]         
         
        #get the significant condidates
        thres_significant = 0.2
        indsSignificant = [i for i, m in enumerate(meritsModified) if m/meritsModified[iZrefExtremum]>thres_significant]
        indBestMerit = 0  
        
        #build the common (all candidates for this source) template vector
        _data = np.zeros((N))
        _data[idxSpcName] = resList.list[k].name.split("_")[1]
        _data[idxMeritStd] = chi2Std
        _data[idxNCandSignificant] = len(indsSignificant)

        #build the data
        for i, kInd in enumerate(indsSignificant):
            data = np.copy(_data) 
            if kInd == iZrefExtremum:
                data[idxFlag] = 1
            else:
                data[idxFlag] = 0
                
            data[idxRedshiftCandidate] = redshifts[kInd]
            data[idxRelRedshiftErrBestCandidate] = (redshifts[kInd]-redshifts[indBestMerit])/(1+redshifts[indBestMerit])
            
            data[idxMeritCandidate] = meritsModified[kInd]
            data[idxRelMeritBestCandidate] = meritsModified[kInd]-meritsModified[indBestMerit]

            #data[idxContIndexColorLya] = chi2.amazed_continuumIndexes[kInd]['color'][0]
            #data[idxContIndexBreakLya] = chi2.amazed_continuumIndexes[kInd]['break'][0]
            data[idxContIndexColorOII] = chi2.amazed_continuumIndexes[kInd]['color'][1]
            data[idxContIndexBreakOII] = chi2.amazed_continuumIndexes[kInd]['break'][1]
            data[idxContIndexColorOIII] = chi2.amazed_continuumIndexes[kInd]['color'][2]
            data[idxContIndexBreakOIII] = chi2.amazed_continuumIndexes[kInd]['break'][2]
            data[idxContIndexColorHa] = chi2.amazed_continuumIndexes[kInd]['color'][3]
            data[idxContIndexBreakHa] = chi2.amazed_continuumIndexes[kInd]['break'][3]
                        
            data_list.append(data)
        
        print("data = {}".format(data_list))
    
    enableExport = 1 
    if enableExport:
        outdir = os.path.join(resList.analysisoutputdir, "svm")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        outFile = os.path.join(outdir, 'svm_data.txt'.format())
        np.savetxt(outFile, data_list)

    
      
def exportContinuumrelevance(resDir, diffthres, spcName=""):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName, methodName="", zrefmin=-1, zrefmax=50)
    if resList.n <1:
        print('No results loaded...')
        return

    resultParser = rp.ResParser(resDir) 
    idxSpcName = 0
    idxStdSpc = 1
    idxStdContinuum = 2

    
    N = 3
    data_list_str = "" #"#name\tstdcont_over_stdspc\n"
    data_list = []
    
    #nres = 100
    nres = resList.n
    for k in range(nres):
        print("\nprocessing result #{}/{}".format(k+1, resList.n))
        
        _spcName = resList.list[k].name
        contRelevancepath = resultParser.getContinuumRelevancePath(_spcName) 
        f = open(contRelevancepath, 'r')
        for l in f:
            if not (l.startswith('#') or 'std_spc' in l):
                print("line is {}".format(l))
                llist = l.split("\t")
                _stdSpc = float(llist[0])
                _stdContinuum = float(llist[1])
                break
        f.close()
                
        print('_stdSpc = {}'.format(_stdSpc))
        print('_stdContinuum = {}'.format(_stdContinuum))
        
        #build the common (all candidates for this source) template vector
#        _data = np.zeros((N))
#        _data[idxSpcName] = _spcName#resList.list[k].name.split("_")[1]
#        _data[idxStdSpc] = _stdSpc
#        _data[idxStdContinuum] = _stdContinuum

        dataStr = "{}\t{}\n".format(_spcName, _stdContinuum)
        data_list_str+=dataStr
        print("data = {}".format(dataStr))
        data_list.append([_stdSpc, _stdContinuum])
    
    enablePlot = 1
    if enablePlot:
        #print("data_list = {}".format(data_list))

###
        histRawData = np.array([a[1]/a[0] for i,a in enumerate(data_list)])
        plt.plot(histRawData, 'x')
        plt.yscale('log')
        plt.grid()
        
        print("histRawData = {}".format(histRawData))
        mini = np.nanmin(histRawData)
        maxi = np.nanmax(histRawData)
        
        print("min={}, max={}".format(mini, maxi))
        nbins = 20
        vectBins = np.linspace(mini, maxi, nbins, endpoint=True)

        mybins = vectBins
                
        print("mybins={}".format(mybins))
        centerbins = [(mybins[k]+mybins[k+1])/2.0 for k in range(len(mybins)-1)]
        widthbins = [(mybins[k+1]-mybins[k]) for k in range(len(mybins)-1)]
        ybins, bin_edges = np.histogram(histRawData, bins=mybins)
     
        fig = plt.figure( "stats", figsize=(15,11))
        if 1:
            width = 0.7 * widthbins[0] #assuming all the bins have the same size
            center = (bin_edges[:-1] + bin_edges[1:]) / 2
            barlist = plt.bar(center, ybins, align='center', width=width)
#            for k in range(len(barlist)):
#                if centerbins[k]<-thres:
#                    barlist[k].set_color('g')   
#                elif centerbins[k]>thres:
#                    barlist[k].set_color('r') 
            
        else:
            plt.plot(centerbins, ybins, 'x-')
        
        nbins = len(ybins)
        plt.xlim(min(mybins), max(mybins))
        #plt.semilogx()
        #plt.xlim([-2, 2])
        #plt.ylim([0, 100])            

        ##bar
        #ind = np.arange(len(OY))
        #plt.plot(xvect, yvect, 'x')
        #ax.set_xscale('log')
        plt.grid(True) # Affiche la grille
        #plt.legend(('cos','sin'), 'uplter right', shadow = True)
        plt.ylabel('Count')
        plt.xlabel('StdContinuum/StdSpectrum'.format())
        plt.show() 
###
    
    
    
    enableExport = 1 
    if enableExport:
        outdir = os.path.join(resList.analysisoutputdir, "continuum")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
            
        outFile = os.path.join(outdir, 'continuumRelevance_data.txt'.format())
        f = open(outFile, "w")
        f.write(data_list_str)
        f.close()

    
    
        
    
def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./resultstat.py """
    parser = optparse.OptionParser(usage=usage)
    #parser.add_option(u"-x", u"--xaxis", help="force the use of indexes for the xaxis (--xaxis=""index""), or enable the use of calibrated wavelength (--xaxis=""angstrom"")",  dest="xAxisType", default="angstrom")
    parser.add_option(u"-d", u"--dir", help="path to the amazed results directory (/output/)",  dest="resDir", default="./output")
    parser.add_option(u"-f", u"--failurediffthres", help="diff threshold for the extraction of the failure spectra",  dest="diffthres", default=-1)


    #parser.add_option(u"-f", u"--diff", help="path to the diff file",  dest="diffFile", default="diff.txt")
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        choiceStr = raw_input("\n\nPlease choose the processing:\n\
        1. Export Full Analysis\n\
        2. Plot Template Missing rate\n\
        3. Plot Closest Z Candidate\n\
        4. Plot Relative Position - Closest Z candidate\n\
        5. Plot Relative Position - 2nd extrema z candidate\n\
        6. Export Best Redshift - with zrange(per template)\n\
        7. Export Sources in z-mag-sfr Bin\n\
        \n\
        10. Export Line Detection Stats (for linematching)\n\
        11. Compare Peak Detection with reference (for linematching)\n\
        12. Compare Line Detection with reference (for linematching)\n\
        13. Export Line Matching Stats (for linematching)\n\
        \n\
        14. Export N Line Stats external.csv file (for lienmodel)\n\
        \n\
        20. Compare failures\n\
        \n\
        30. Plot 2D Combination Merit Coeff map\n\
        31. Estimate amazed03 coeff map\n\
        32. Export amazed03 SVM table\n\
        35. Plot continuum indexes (only for linemodel)\n\
        36. Export SVM table\n\
        37. Export continuum Relevance\n\
        \n\
        40. Plot precision histogram\n\
        41. Plot velocity error histogram\n\
        \n")
        choice = int(choiceStr)
        
        #non implemented yet
        #plotChi2XYPlane(options.resDir)
        
        #plotReducedZcandidates(options.resDir)
            
            
        if choice == 1:
            exportAnalysis(options.resDir)
        elif choice == 2:
            plotTplMissingRate(options.resDir)
        elif choice == 3:
            extremaTypeStr = raw_input("Please enter the extrema type : choices = raw, continuum, nocontinuum, corr, linemodel, linemodeltplshape :")
            if not (extremaTypeStr == "raw" or extremaTypeStr == "continuum" or extremaTypeStr == "nocontinuum" or extremaTypeStr == "corr" or extremaTypeStr == "linemodel" or extremaTypeStr == "linemodeltplshape"):
                print("extrema type not successfully, aborting")
                return
                
            extrChoiceStr = raw_input("\n\nPlease enter the number of extrema to be considered...\n")
            nextrInput = int(extrChoiceStr)
            plotClosestZCandidates(options.resDir, float(options.diffthres), nextrInput, extremaTypeStr)
        elif choice == 4:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
            extremaTypeStr = raw_input("Please enter the extrema type : choices = raw, continuum, nocontinuum, corr, linemodel, linemodeltplshape :")
            if not (extremaTypeStr == "raw" or extremaTypeStr == "continuum" or extremaTypeStr == "nocontinuum" or extremaTypeStr == "corr" or extremaTypeStr == "linemodel" or extremaTypeStr == "linemodeltplshape"):
                print("extrema type not successfully, aborting")
                return
            extrChoiceStr = raw_input("\n\nPlease enter the number of extrema to be considered...\n")
            nextrInput = int(extrChoiceStr)
            
            plotRelativePosClosestZref(options.resDir, spcName, nextrInput, extremaTypeStr)
        elif choice == 5:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
            plotRelativePosSecondBestExtrema(options.resDir, float(options.diffthres), spcName)
        
        elif choice == 6:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
            extremaTypeStr = raw_input("Please enter the extrema type : choices = raw, nocontinuum, continuum, corr, linemodel, linemodeltplshape :")
            if not (extremaTypeStr == "raw" or extremaTypeStr == "continuum" or extremaTypeStr == "nocontinuum" or extremaTypeStr == "corr" or extremaTypeStr == "linemodel" or extremaTypeStr == "linemodeltplshape"):
                print("extrema type not successfully, aborting")
                return
            enableZrangeFilterStr = raw_input("Do you want to use the hardcoded zrange values to filter each z candidates ? (y, n) :")
            if not (enableZrangeFilterStr == "y" or enableZrangeFilterStr == "n"):
                print("enableZrangeFilterStr not successfully retrieved from user input, aborting")
                return
            else:
                if enableZrangeFilterStr == "n":
                    enableZrangeFilter=False
                if enableZrangeFilterStr == "y":
                    enableZrangeFilter=True
                    
            exportBestRedshiftWithZRangePerTemplate(options.resDir, float(options.diffthres), chi2Type=extremaTypeStr, spcName=spcName, enableZrangeFilter=enableZrangeFilter)
                
        elif choice == 7:
            zrefmin = 4.5
            zrefmax = 5.0
            print("using zrange = {:<10}{:<10}".format(zrefmin, zrefmax))
            magrefmin = 21.0
            magrefmax = 22.0
            print("using magrange = {:<10}{:<10}".format(magrefmin, magrefmax))
            sfrrefmin = 10.0
            sfrrefmax = 100.0
            print("using sfrrnage = {:<10}{:<10}".format(sfrrefmin, sfrrefmax))
            printSourcesInZMagSfrBin(options.resDir, float(options.diffthres), zrefmin=zrefmin, zrefmax=zrefmax, magrefmin=magrefmin, magrefmax=magrefmax, sfrrefmin=sfrrefmin, sfrrefmax=sfrrefmax)
            
        elif choice == 10:
            exportLineDetectionStats(options.resDir)
        elif choice == 11:
            refresDir = raw_input("Enter a reference dir. to compare the peak detection results ? (press enter to skip) :")
            refresDir = refresDir.strip()
            refresDir = refresDir.replace("'", "")
            if not (refresDir == ""):
                comparePeakDetection(options.resDir, refresDir)
        elif choice == 12:
            refresDir = raw_input("Enter a reference dir. to compare the detection results ? (press enter to skip) :")
            refresDir = refresDir.strip()
            refresDir = refresDir.replace("'", "")
            if not (refresDir == ""):
                compareLineDetection(options.resDir, refresDir)
        elif choice == 13:
            exportLineMatchingStats(options.resDir)
                            
        elif choice == 14:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
           
            exportNLinesExternal(options.resDir, float(options.diffthres), spcName)
            
        
        elif choice == 20:
            refresDir = raw_input("Enter a reference dir. to compare the failures ? (press enter to skip) :")
            refresDir = refresDir.strip()
            refresDir = refresDir.replace("'", "")
            if refresDir == "":
                print("ERROR: empty reference result directory: aborting...")    
            else:
                diffthreshold = 0.01
                zrefmin = -1
                zrefmax = 50.0
                print("INFO: using default diffthreshold={}, and zrange=[{} {}]".format(diffthreshold, zrefmin, zrefmax)) 
                WarningKeyStr = raw_input("Press any key to continue...".format())
        
                compareFailures(options.resDir, refresDir, diffthreshold=diffthreshold, zrefmin=zrefmin, zrefmax=zrefmax) 
                
        elif choice == 30:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
            methodName = ""
            methodStr = raw_input("Do you want to enter a method name to filter the results ? (press enter to skip) :")
            if not (methodStr == "No" or methodStr == "no"):
                methodName = methodStr
            plotChi2CombinationCoeff2DMap(options.resDir, float(options.diffthres), spcName, methodName)
                            
        elif choice == 31:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr

            estimateCombinationCoeffMap(options.resDir, float(options.diffthres), spcName)
                                        
        elif choice == 32:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr

            exportAmazed03SVMTable(options.resDir, float(options.diffthres), spcName)
            
                
        elif choice == 35:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
           
            plotContinuumIndexes(options.resDir, float(options.diffthres), spcName)
            

        elif choice == 36:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
           
            exportSVMTable(options.resDir, float(options.diffthres), spcName)
            
        elif choice == 37:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
           
            exportContinuumrelevance(options.resDir, float(options.diffthres), spcName)
            
        elif choice == 40:           
            plotPrecisionHist(options.resDir)            
        elif choice == 41:           
            plotVelocityError(options.resDir)
            
        else:    
            print("Error: invalid entry, aborting...")
    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()
   
if __name__ == '__main__':
    Main( sys.argv )

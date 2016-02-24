# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 10:49:18 2015

@author: aschmitt
"""
import sys
import os
import inspect
import optparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import math
import time

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

class Result(object):
    def __init__(self, name, zref, zcalc, zdiff, chi2, chi2nc, chi2PerTplZcalc, chi2PerTplZref, chi2ncPerTplZref):
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
    def __init__(self, dir, diffthreshold=-1, opt='full', spcName="", methodName=""):
        self.logTagStr = "ResultList"
        self.dir = dir
        self.name = os.path.basename(self.dir)
        self.diffthreshold = diffthreshold  
        self.opt = opt
        self.statsdir = ""
        self.list = [] 
        self.n = -1;

        self.analysisoutputdir = ""
        
        
        self.load(spcName=spcName, methodName=methodName)
    
    def load(self, spcName="", methodName=""):
        s = rp.ResParser(self.dir)
        print(s) 
        self.statsdir = s.getStatsDirPath()
        self.analysisoutputdir = os.path.join(self.statsdir, "analysis")
        diffthresStr = "diffthres_{}".format(self.diffthreshold)
        if(self.diffthreshold<0):
            diffthresStr = "diffthres_none"
        self.analysisoutputdir = os.path.join(self.analysisoutputdir, diffthresStr)
        if not os.path.exists(self.analysisoutputdir):
            os.makedirs(self.analysisoutputdir)
    
        n = s.getDiffSize() 
        #n = 4
        print('Number of spectra is: {0}'.format(n))
    
        self.n = 0;
        for x in range(0, n):
            line = s.getDiffLine(x)
            name = line[0] 
            zref = line[2]
            zcalc = line[4]
            zdiff = line[4]-line[2]
            #print('zcalc is: {0}'.format(zcalc)) 
            method = line[6]
            
            if 1:
                print('name is: {0}'.format(name))
                print('zref is: {0}'.format(zref))
                print('zcalc is: {0}'.format(zcalc)) 
                print('self.diffthreshold is: {0}'.format(self.diffthreshold))
                print('method is: {0}'.format(method)) 
            if self.opt=='full':
                chi2 = s.getChi2Val(name, zcalc, "", chi2type="raw")
                chi2nc = s.getChi2Val(name, zcalc, "", chi2type="nocontinuum")
                chi2PerTplZcalc = s.getChi2List(name, zcalc, chi2type="raw")
                chi2PerTplZref = s.getChi2List(name, zref, chi2type="raw")
                chi2ncPerTplZref = [];#s.getChi2List(name, zref, chi2type="nocontinuum")
            else:
                chi2 = -1
                chi2nc = -1
                chi2PerTplZcalc = []
                chi2PerTplZref = []
                chi2ncPerTplZref = []
            if (spcName=="" or spcName==name) and (methodName=="" or methodName==method):   
                if abs(zdiff)>=self.diffthreshold:
                    res = Result(name, zref, zcalc, zdiff, chi2, chi2nc, chi2PerTplZcalc, chi2PerTplZref, chi2ncPerTplZref)
                    self.list.append(res)
                    self.n +=1
            #print('{0}/{1} results loaded: last chi2 = {2}'.format(self.n,x, chi2))
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
    
    def getZCandidatesFromAmazedChi2Extrema(self, indice=0, chi2Type = "raw", nextrema = 5):
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
                print("Problem while retrieving chi2 filepath..")
            else:
                print("using Chi2 file path : ".format(filepath))
            chi2 = chisq.ResultChisquare(filepath)
            
            redshiftstpl = chi2.amazed_extrema
            if len(redshiftstpl)<1:
                return [-1],[-1]
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
            redshifts = [redshiftsUnsorted[b] for b in sortId]
            merits = [meritsUnsorted[b] for b in sortId]
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
                redshifts, merits = self.getZCandidatesFromAmazedChi2Extrema(x, chi2Type, nextrema)
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
                    zabsdiff[x] = abs(self.list[x].zref-z)
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
                
                fout.write( "idx" + '\t' + "spc" + '\t' + 'zclosest' + '\t' + "rank" + '\t' + "abszdiff"  + '\n')
                foutFailures.write( "idx" + '\t' + "spc" + '\t' + 'zclosest' + '\t' + "rank" + '\t' + "abszdiff"  + '\n')
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
            plt.xlabel('abs(zref - closest z candidate)')
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
        
    def getSecondExtremaZcandidateRelativePositionList(self, chi2Type="raw", extremaType="amazed", enablePlot=True, enableExport=False):
        zrelativemerit = []
        diffvect = []
        n = self.n
        #n = 10#self.n
        nextrema = 15
                
        
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

    
    def getChi2LinCombinationCoeff2DMap(self, indice=0, enablePlot=0, chi2dontloadThres=-1, zthres=0.01):       
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
#        for i in range(len(chi2_continuum)):
#            chi2_continuum[i] = chi2_continuum[i]#/nspcSamples
#            chi2_continuum[i] = chi2_continuum[i]*nspcSamples
#            if max_chi2_all<chi2_continuum[i] and chi2_continuum[i]<1e20:
#                max_chi2_all = chi2_continuum[i]
        chi2type = "nocontinuum"
        z_nc, chi2_nc = s.getChi2MinValList(spcName, chi2type, dontloadThres=chi2dontloadThres) 
#        for i in range(len(chi2_nc)):
#            chi2_nc[i] = chi2_nc[i]#/nspcSamples
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
            ax.plot(z_nc, chi2_nc, label='no continuum')
            ax.plot(z_lm, chi2_lm, label='linemodel', linestyle = "dashed")
            ax.plot(z_continuum, chi2_continuum, label='continuum')
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
        coeff_nc = 100
        coeff_continuum_min = 0.0
        coeff_continuum_max = 99.0
        coeff_continuum_step = 5.0
        n_continuum = int((coeff_continuum_max-coeff_continuum_min+1)/float(coeff_continuum_step))
        
        coeff_lm_min = 0.0
        coeff_lm_max = 99.0
        coeff_lm_step = 5.0
        n_lm = int((coeff_lm_max-coeff_lm_min+1)/float(coeff_lm_step))
        merit = np.zeros((nz))
        _map = np.zeros((n_continuum, n_lm))
        
        #print("n_continuum = {}".format(n_continuum))
        #print("n_lm = {}".format(n_lm))
        for icontinuum in range(n_continuum):
            print "."
            coeff_continuum = coeff_continuum_min+icontinuum*coeff_continuum_step
            for ilinemodel in range(n_lm):
                coeff_lm = coeff_lm_min+ilinemodel*coeff_lm_step
                #print("this solution uses : coeff_continuum={}, coeff_lm={}".format(coeff_continuum, coeff_lm) )
                for i in range(nz):
                    merit[i] = coeff_nc*chi2_nc[i] + coeff_continuum*chi2_continuum[i] + coeff_lm*chi2_lm[i]
                izbest = np.argmin(merit)
                zbest = z_nc[izbest]
                if(np.abs(zbest - self.list[indice].zref ) < zthres):
                    _map[icontinuum,ilinemodel] = 1
                    #print("this solution works : z={}, zref={}".format(zbest, self.list[indice].zref) )
                else:          
                    _map[icontinuum,ilinemodel] = 0
                    #print("this solution fails : z={}, zref={}".format(zbest, self.list[indice].zref) )
        return _map
    
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

        removeStrRef = "_F"
        #removeStrRef = "_FILT_MGv0_c"
        #removeStrRef = "_FILT3dec2015"
        #removeStrThis = "_FILT_MGv1_c"
        removeStrThis = "_F"
        #removeStr2 = "_FILT"        
        
        n = self.n
        nref = refreslist.n
        for x in range(0,n):
            #print("\n")
            #print("Spc {}/{}".format(x+1, self.n))
            spcName = self.list[x].name
            spcName = spcName.replace(removeStrThis, "")
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
                
        for xref in range(0,nref):
            #print("\n")
            #print("Spc {}/{}".format(x+1, self.n))
            spcNameRef = refreslist.list[xref].name
            spcNameRef = spcNameRef.replace(removeStrRef, "")
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
                
        print("\n\n")
        print("########## Failures Comparison: ##########")
        print("N same failures = {}".format(len(failures_both)))
        print("\n")
        print("N total failures this = {}".format(n))
        print("N failures only this = {}".format(len(failures_onlyThis)))
        for x in range(0,len(failures_onlyThis)):
            print("\t{}".format(failures_onlyThis[x]))
        print("\n")
        print("N total failures ref = {}".format(nref))
        print("N failures only ref = {}".format(nref-len(failures_both)))
        for x in range(0,len(failures_onlyRef)):
            print("\t{}".format(failures_onlyRef[x]))
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
                thisStr = thisStr +"\t{}\n".format(failures_onlyThis[x])
            f.write(thisStr+"\n")  
            refStr = "Failures for REF resultset = {} :\n".format(refreslist.name)
            refStr = refStr + "N total failures REF resultset = {}\n".format(nref)
            refStr = refStr + "N failures only REF resultset = {}\n".format(nref-len(failures_both))        
            for x in range(0,len(failures_onlyRef)):
                refStr = refStr +"\t{}\n".format(failures_onlyRef[x])
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
    
#    relposraw = resList.getClosestZcandidateRelativePositionList("raw", "amazed", enablePlot=True, enableExport=False, nextrema = nextrema)
#    zrelativemeritList.append(relposraw)
#    relposnocont = resList.getClosestZcandidateRelativePositionList("nocontinuum", "amazed", enablePlot=True, enableExport=False, nextrema = nextrema)
#    zrelativemeritList.append(relposnocont)
#    relposcorr = resList.getClosestZcandidateRelativePositionList("corr", "amazed", enablePlot=True, enableExport=False, nextrema = nextrema)
#    zrelativemeritList.append(relposcorr)
    relposlinemodel = resList.getSecondExtremaZcandidateRelativePositionList("linemodel", "amazed", enablePlot=True, enableExport=True)
    zrelativemeritList.append(relposlinemodel)
    
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
      
def plotChi2LinCombinationCoeff2DMap(resDir, diffthres, spcName="", methodName="", enableExport=True):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=diffthres, opt='brief', spcName=spcName, methodName=methodName)
    if resList.n <1:
        print('No results loaded...')
        return
    cumulcoeffmap = np.zeros((20,20))
    
    ## parameters :
    zthres = 0.01
    if 0: 
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
        coeffmap = resList.getChi2LinCombinationCoeff2DMap(k, enablePlot=enablePlot, chi2dontloadThres=chi2dontloadThres,  zthres=zthres)
        print("coeffmap = {}".format(coeffmap))
        cumulcoeffmap = np.add(cumulcoeffmap, coeffmap)
        print("cumulcoeffmap = {}".format(cumulcoeffmap))
        #print("coeffmap shape = {}".format(coeffmap.shape))  
        plt.clf()
        plt.close()

        fig = plt.figure('dtreeb lincomb coeff map', figsize=(9, 8))
        ax = fig.add_subplot(111)
        cmap = plt.get_cmap('RdYlGn') 
        ncolors = 20 #min(20,k+2)
        cmap = cmap_discretize(cmap, ncolors) 
        
        i = ax.matshow(np.transpose(cumulcoeffmap), interpolation='nearest', aspect='equal', cmap=cmap)
        plt.xlabel('continuum coeff')
        plt.ylabel('lm coeff')
        name1 = "linear combination coeff map \nzThreshold = {}\n(processed idx={}/{})".format(zthres, k+1, resList.n)
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
            outdir = os.path.join(resList.analysisoutputdir, "dtreeb_lincomb_coeffmap_method{}".format(methodName))
            if not os.path.exists(outdir):
                print("creating outputdir {}".format(outdir))
                os.makedirs(outdir)   
            outFigFile = os.path.join(outdir, 'cumulcoeffmap.png')
            plt.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
            outtxtFile = os.path.join(outdir, 'cumulcoeffmap.txt')
            np.savetxt(outtxtFile, cumulcoeffmap)
   
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
   
def compareFailures(resDir, refresdir):
    print('using amazed results full path: {0}'.format(resDir))
    resList = ResultList(resDir, diffthreshold=0.01, opt='brief')
    refresList = ResultList(refresdir, diffthreshold=0.01, opt='brief')
    resList.getFailuresComparison(refresList)
  

    
    
    
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
        3. Plot Closest Z Candidats\n\
        4. Plot Relative Position - Closest Z candidate\n\
        5. Plot Relative Position - 2nd extrema z candidate\n\
        \n\
        10. Export Line Detection Stats\n\
        11. Compare Peak Detection with reference\n\
        12. Compare Line Detection with reference\n\
        13. Export Line Matching Stats\n\
        \n\
        20. Compare failures\n\
        \n\
        30. Plot 2D Linear Comb. Merit Coeff map\n\
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
            extremaTypeStr = raw_input("Please enter the extrema type : choices = raw, nocontinuum, corr, linemodel :")
            if not (extremaTypeStr == "raw" or extremaTypeStr == "nocontinuum" or extremaTypeStr == "corr" or extremaTypeStr == "linemodel"):
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
            extremaTypeStr = raw_input("Please enter the extrema type : choices = raw, nocontinuum, corr, linemodel :")
            if not (extremaTypeStr == "raw" or extremaTypeStr == "nocontinuum" or extremaTypeStr == "corr" or extremaTypeStr == "linemodel"):
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
        
        elif choice == 20:
            refresDir = raw_input("Enter a reference dir. to compare the failures ? (press enter to skip) :")
            refresDir = refresDir.strip()
            refresDir = refresDir.replace("'", "")
            if not (refresDir == ""):
                compareFailures(options.resDir, refresDir) 
                
        elif choice == 30:
            spcName = ""
            spcStr = raw_input("Do you want to enter a spectrum name to filter the results ? (press enter to skip) :")
            if not (spcStr == "No" or spcStr == "no"):
                spcName = spcStr
            methodName = ""
            methodStr = raw_input("Do you want to enter a method name to filter the results ? (press enter to skip) :")
            if not (methodStr == "No" or methodStr == "no"):
                methodName = methodStr
            plotChi2LinCombinationCoeff2DMap(options.resDir, float(options.diffthres), spcName, methodName)
            
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

# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 11:44:49 2015

@author: aschmitt
"""
import os, sys
import math
import optparse

import matplotlib as mpl
mpl.use('Qt5Agg')

import matplotlib.pyplot as pp
#pp.style.use('ggplot')
#import seaborn as sns
#sns.set_context("poster")
#sns.set_style("whitegrid")
import scipy.interpolate
import numpy as np

import catalog as ctlg

class ResultChisquare(object):
    def __init__(self, spath, stype='undef', dontloadThres=1e31, label=""):
        self.logTagStr = "ResultChisquare"
        self.spath = spath
        tplpath = os.path.dirname(spath)
        tplname = os.path.basename(tplpath)      
        spcpath = os.path.dirname(tplpath)
        spcname = os.path.basename(spcpath)  
        self.name = spcname + "\n" + tplname + "\n" + os.path.basename(spath)
        if label=="":
            self.label = self.name
        else:
            self.label = label
        self.stype = stype
        self.n = -1
        self.xvect = -1
        self.yvect = -1
        self.ysum = 0
        self.forcePlotXIndex = False
        self.amazed_extrema = []
        self.amazed_logarea = []
        self.amazed_sigmaz = [] 
        self.amazed_fitamplitude = []
        self.amazed_continuumIndexes = [] 
        self.amazed_StrongELSNR = []     
        
        #self.cpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/RayCatalogs/raycatalogamazedvacuum.txt"
        self.cpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/linecatalogs/linecatalogamazedair_B.txt"
        #self.cpath = "/home/aschmitt/data/pfs/pfs_reallyjustline/amazed/RayCatalogs/raycatalogamazedair_custompfs2.txt"
        self.load(dontloadThres=dontloadThres)
        
    def load(self, dontloadThres=1e31):
        filename = self.spath
        wave = []
        flux = []
        amazed_continuumIndexesColor = []
        amazed_continuumIndexesBreak = []
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                data = lineStr.split("\t")
                #data = lineStr.split(" ")
                #data = [r for r in data if r != '']
                #print len(data)
                if(len(data) >=2):               
                    if(float(data[1])<dontloadThres or dontloadThres<=0):
                        # fill the list                
                        #print data[0]
                        #print data[1]
                        wave.append(float(data[0]))
                        flux.append(float(data[1]))
            elif not lineStr.find("Extrema for z = {") == -1:
                beg = lineStr.find("{")
                end = lineStr.find("}")  
                dataStr = lineStr[beg+1:end]
                data = dataStr.split("\t")
                data = [d for d in data if len(d.strip(" "))>0]
                for d in data:
                    self.amazed_extrema.append(float(d))
            elif not lineStr.find("SigmaZ for each extrema = {") == -1:
                beg = lineStr.find("{")
                end = lineStr.find("}")  
                dataStr = lineStr[beg+1:end]
                data = dataStr.split("\t")
                data = [d for d in data if len(d.strip(" "))>0]
                for d in data:
                    self.amazed_sigmaz.append(float(d))
            elif not lineStr.find("LogArea for each extrema = {") == -1:
                beg = lineStr.find("{")
                end = lineStr.find("}")  
                dataStr = lineStr[beg+1:end]
                data = dataStr.split("\t")
                data = [d for d in data if len(d.strip(" "))>0]
                for d in data:
                    self.amazed_logarea.append(float(d))
            elif not lineStr.find("Extrema FitAmplitudes = {") == -1:
                beg = lineStr.find("{")
                end = lineStr.find("}")  
                dataStr = lineStr[beg+1:end]
                data = dataStr.split("\t")
                data = [d for d in data if len(d.strip(" "))>0]
                for d in data:
                    self.amazed_fitamplitude.append(float(d))
            elif not lineStr.find("ContinuumIndexes Color for each extrema = {") == -1:
                beg = lineStr.find("{")
                end = lineStr.find("}")  
                dataStr = lineStr[beg+1:end]
                strExtrema = dataStr.split(">")
                strExtrema = [d.replace("<", "") for d in strExtrema]
                #print("strExtrema = {}".format(strExtrema))
                dataExtrema = [d for d in strExtrema if len(d.strip(" "))>0]
                for de in dataExtrema:
                    #print("de = {}".format(de))
                    extremaContIndexes = []
                    data = de.split("\t")
                    data = [d for d in data if len(d.strip(" "))>0]
                    for d in data:
                        #print("d = {}".format(d))
                        extremaContIndexes.append(float(d))
                    amazed_continuumIndexesColor.append(extremaContIndexes)
            elif not lineStr.find("ContinuumIndexes Break for each extrema = {") == -1:
                beg = lineStr.find("{")
                end = lineStr.find("}")  
                dataStr = lineStr[beg+1:end]
                strExtrema = dataStr.split(">")
                strExtrema = [d.replace("<", "") for d in strExtrema]
                #print("strExtrema = {}".format(strExtrema))
                dataExtrema = [d for d in strExtrema if len(d.strip(" "))>0]
                for de in dataExtrema:
                    #print("de = {}".format(de))
                    extremaContIndexes = []
                    data = de.split("\t")
                    data = [d for d in data if len(d.strip(" "))>0]
                    for d in data:
                        #print("d = {}".format(d))
                        extremaContIndexes.append(float(d))
                    amazed_continuumIndexesBreak.append(extremaContIndexes)
            elif not lineStr.find("StrongELSNR for each extrema =") == -1:
                beg = lineStr.find("{")
                end = lineStr.find("}")  
                dataStr = lineStr[beg+1:end]
                data = dataStr.split("\t")
                data = [d for d in data if len(d.strip(" "))>0]
                for d in data:
                    self.amazed_StrongELSNR.append(float(d))
                
        f.close()
        self.n = len(wave)
        #print('len wave = {0}'.format(self.n))
        #---- default xaxis index array
        self.xvect = list(range(0,self.n))
        self.yvect = list(range(0,self.n)) 
        for x in range(0,self.n):
            self.xvect[x] = wave[x]
            #self.yvect[x] = math.exp(-flux[x]/2.0)
            self.yvect[x] = flux[x]
            
        nconti = len(amazed_continuumIndexesColor)
        #print("amazed_continuumIndexesColor size = {}".format(len(amazed_continuumIndexesColor)))
        #print("amazed_continuumIndexesBreak size = {}".format(len(amazed_continuumIndexesBreak)))
        for x in range(0,nconti):
            self.amazed_continuumIndexes.append({'color': amazed_continuumIndexesColor[x], 'break': amazed_continuumIndexesBreak[x]})
            
            

        #print("Extrema found : {0}".format(self.amazed_extrema))
        print("loaded chi2 : {0}".format((filename)))
        
        
        
    def __str__(self):
        a = "\nChisquare: {0}\n".format(self.name)
        a = a + ("    type = {0}\n".format(self.stype))
        a = a + ("    n = {0}\n".format(self.n))
        a = a + ("    flux min = {0}\n".format(self.getFluxMin()))
        a = a + ("    flux median = {0}\n".format(self.getFluxMedian()))
        a = a + ("\n")
        
        for z in range(len(self.amazed_extrema)):
            extrema = -1
            try:
                extrema = self.amazed_extrema[z]
            except:
                pass
            logarea = -1
            try:
                logarea = self.amazed_logarea[z]
            except:
                pass
            sigmaz = -1
            try:
                sigmaz = self.amazed_sigmaz[z]
            except:
                pass
            strongELSNR = -1
            try:
                strongELSNR = self.amazed_StrongELSNR[z]
            except:
                pass
            
            a = a + ("    zcandidate_{}\tz = {}\t\tlogarea = {}\t\tsigmaz = {}\t\tstrongELSNR = {}\n".format( z, extrema, logarea, sigmaz, strongELSNR )) 
            
        a = a + ("\n")
        
        a = a + ("    contIndexes = {:<20}\t{:<20}\t{:<20}\t{:<20}\n".format("Lya", "OII", "OIII", "Halpha"))
        for z in range(len(self.amazed_continuumIndexes)):
            dataColor = self.amazed_continuumIndexes[z]['color']
            dataBreak = self.amazed_continuumIndexes[z]['break']
            a = a + ("    ci_color_{} = ".format(z))
            for d in dataColor:
                a = a + ("{:<20}\t".format(d ))
            a = a + ("\n") 
            a = a + ("    ci_break_{} = ".format(z))
            for d in dataBreak:
                a = a + ("{:<20}\t".format(d ))
            a = a + ("\n")   
        
        
        return a
        
    def getXYSortedByY(self):
        x = self.xvect
        y = self.yvect
        
        sortId=np.argsort(y)
        
        x=[x[s] for s in sortId]
        y=[y[s] for s in sortId]
        return [x, y]
      
    def getExtrema(self, idx):
        if idx < len(self.amazed_extrema):
            return self.amazed_extrema[idx]
        else:
            return -1.0;      
    
    def plot(self, showContinuumEstimate=False, showExtrema=False, showAmbiguities=False, enablePlot=True, exportPath="", enableReturnFig=False):
        showStrongElPriors = True
        thresShortMeritCurve = 25
        #find limits
        cmin = +1e6;
        cmax = -1e6;
        thres = 1e6
        for x in range(self.n):
            if cmin >  self.yvect[x] and self.yvect[x]>-thres:
                cmin = self.yvect[x]
            if cmax <  self.yvect[x] and self.yvect[x]<thres:
                cmax = self.yvect[x]
        #overrride cmax
        if self.n > thresShortMeritCurve:     
            #cmax = self.getMeanValue()
            cmax = self.getFluxMedian()
            
            
        range_ori = (cmax-cmin)
        cmax = cmax+0.05*range_ori
        cmin = cmin-0.2*range_ori
        
        
        fig = pp.figure("chi2")
        titleStr = self.name
        if self.n < thresShortMeritCurve: #enable marker when there are only few samples
            marker='o'
            linestyle=''
        else:
            marker=""
            linestyle='-'
        if not self.forcePlotXIndex:
            pp.plot(self.xvect, self.yvect, linestyle=linestyle, marker=marker)
        else:
            pp.plot(self.yvect, linestyle=linestyle, marker=marker)

        if showContinuumEstimate:
            chi1cont = self.getSplineContinuum()
            if not self.forcePlotXIndex:
                pp.plot(self.xvect, chi1cont, 'r')
            else:
                pp.plot(chi1cont, 'r') 
        if showExtrema:
            linesx = self.amazed_extrema
            print("chi2: loaded extrema = {}".format(linesx))
            for k in range(len(linesx)):
                x = linesx[k]
                minivect = [np.abs(a-x) for a in self.xvect]
                if len(minivect)>0:
                    idx = np.argmin(minivect)
                    chimin = self.yvect[idx]*1.0
                    chimax = self.yvect[idx]*1.0
                    pp.plot((x, x), (chimin, chimax) , 'ro', label="Extremum_{}".format(k) ) 
            if showStrongElPriors:
                linesy = self.get_StrongELSnrPrior()
                print("chi2: loaded priors = {}".format(linesy))
                if len(linesy) == len(linesx):
                    for k in range(len(linesx)):
                        x = linesx[k]
                        minivect = [np.abs(a-x) for a in self.xvect]
                        if len(minivect)>0:
                            idx = np.argmin(minivect)
                            chimin = self.yvect[idx] - linesy[k]
                            chimax = self.yvect[idx]
                            pp.plot((x, x), (chimin, chimax) , 'g-+', label="SnrPrior_{}".format(k) )
        if showAmbiguities and len(self.amazed_extrema)>0:
            zrefamb = self.amazed_extrema[0]
            #zrefamb = 1.064
            #zrefamb = 0.691323
            titleStr = titleStr + "\nambiguities for z={}".format(zrefamb)
            c = ctlg.Catalog(self.cpath)
            print(c) 
            ambig = c.getAmbiguityZ(zrefamb, self.xvect[0], self.xvect[len(self.xvect) -1], lineTypeFilter=-1, lineForceFilterReference="S", lineForceFilterFail="S")
            linesx = [a["z"] for a in ambig]
            ambig_name = ["{}_{}/{}_{}".format(a["reftype"], a["refname"],a["failtype"], a["failname"]) for a in ambig]
            
            print("chi2: loaded ambig.  n = {}".format(len(linesx)))
            for k in range(len(linesx)):
                x = linesx[k]
                minivect = [np.abs(a-x) for a in self.xvect]
                idx = np.argmin(minivect)
                chimin = self.yvect[idx] - (self.yvect[idx] - cmin)*0.2
                chimax = cmin 
                pp.plot((x, x), (chimin, chimax) , 'k-', label="Extremum_{}".format(k) )
                pp.text(x, chimax, '{0}'.format(ambig_name[k]))

        pp.grid(True) # Affiche la grille
        #pp.legend(('cos','sin'), 'upper right', shadow = True)
        if not self.forcePlotXIndex:
             pp.xlabel('z')
        else:
            pp.xlabel('index')
        pp.ylabel('y')
        pp.title(titleStr) # Titre  
        pp.ylim([cmin,cmax])
        #pp.savefig('ExempleTrace') # sauvegarde du fichier ExempleTrace.png
        
        
        if not enableReturnFig:
            fig.tight_layout()
            if enablePlot:
                pp.show()
            else:
                outFigFile = os.path.join(exportPath, 'chi2_{}.png'.format(self.stype))
                #pp.savefig( outFigFile, bbox_inches='tight')
                pp.savefig( outFigFile)
                pp.close()
                pp.clf()
        else:
            return fig
            
        print('\n')
        
    def plotCompare(self, other_spc):
        #find limits
        cmin = +1e6;
        cmax = -1e6;
        thres = 1e6
        for x in range(self.n):
            if cmin >  self.yvect[x] and self.yvect[x]>-thres:
                cmin = self.yvect[x]
            if cmax <  self.yvect[x] and self.yvect[x]<thres:
                cmax = self.yvect[x]
             
        
        selfColor = 'g'#, '0.5'
        fig, ax1 = pp.subplots()
        lns1 = ax1.plot(other_spc.xvect, other_spc.yvect, color=selfColor,linewidth = 2.0, label=other_spc.label)
        ax1.set_ylabel('Least-Square Error', color=selfColor)
        for tl in ax1.get_yticklabels():
            tl.set_color(selfColor)
    
        otherColor = 'b'
        if 0:
            ax2 = ax1.twinx() 
            if 1:
                ax2.plot(self.xvect, self.yvect, otherColor, linestyle = "-", linewidth = 1.0, label=self.label)
                ax2.set_ylabel('Least-Square Error', color=otherColor)
            else:  
                lns2 = ax2.plot(self.xvect, self.yvect, otherColor, linestyle = "dashed", linewidth = 1.0, label="continuum")# label=self.label)
                
                #chi2path3 = "/home/aschmitt/data/vuds/VUDS_flag3_4/amazed/output/sc_530036900_F53P001_join_seq.fits_60_1_atm_clean/COMBINE-ave-Lya-emstr-AND-StarBurst3.txt/chisquare2solve.chisquare.csv"
                #chi2path3 = "/run/media/aschmitt/al1passport/data_transfert/res_20160511_chisquareRaw/sc_530002417_F53P002_join_A_93_1_atm_clean/NEW-Im-extended.txt/chisquare2solve.chisquare.csv"            
                #chi2path3 = "/run/media/aschmitt/al1passport/data_transfert/res_20160511_chisquareRaw/sc_530004745_F53P002_join_A_46_1_atm_clean/COMBINE-ave-Lya-abs-AND-StarBurst1.txt/chisquare2solve.chisquare.csv"            
                chi2path3 = "/run/media/aschmitt/al1passport/data_transfert/res_20160511_chisquareRaw/sc_530011852_F53P002_join_A_56_1_atm_clean/COMBINE-ave-Lya-emstr-AND-Scd.txt/chisquare2solve.chisquare.csv"            
                s3 = ResultChisquare(chi2path3)
                lns3 = ax2.plot(s3.xvect, s3.yvect, otherColor, linestyle = "-", linewidth = 2.0, label="raw")# label=s3.label)
                ax2.set_ylabel('Least-Square Error', color=otherColor)
            
            
            for tl in ax2.get_yticklabels():
                tl.set_color(otherColor)
        else:
            ax1.plot(self.xvect, self.yvect, otherColor, linestyle = "", marker='v', markersize=9, linewidth = 2.0, label=self.label)

        pp.grid(True) # Affiche la grille
        if 0:
            # added these three lines
            lns = lns1+lns2+lns3
            labs = [l.get_label() for l in lns]
            pp.legend(lns, labs, loc=0)
        if not self.forcePlotXIndex:
            ax1.set_xlabel('z')
        else:
            ax1.set_xlabel('index')
        #pp.ylabel('y')
        #pp.title(self.name) # Titre
        #pp.ylim([cmin,cmax])
        #pp.savefig('ExempleTrace') # sauvegarde du fichier ExempleTrace.png
            
        pp.legend(loc=0)
        pp.show()
        print('\n')

    def getSplineContinuum(self):
        x = np.copy(self.xvect)
        y = np.copy(self.yvect)
        sz = 0.13/(x[1]-x[0])
        sp = scipy.interpolate.UnivariateSpline(x, y, s=sz)
        #pp.plot(x,sp(x))
        #print "splerr", scipy.interpolate.ssqe(sp(x), s, npts)
        #pp.show()
        return sp(x)         
               
    def getFluxMedian(self):
        return np.median(self.yvect)
        
    def getMeanValue(self):
        return np.mean(self.yvect)
        
    def getFluxStd(self):
        return np.std(self.yvect)
        
    def getFluxMin(self):
        return min(self.yvect)
        
    def getFluxMax(self):
        return max(self.yvect)
    
    def getBestZ_byFluxMin(self, zmin, zmax):
        imin = self.getZIndex(zmin)
        imax = self.getZIndex(zmax)
        try:
            iFluxMin = np.argmin(self.yvect[imin:imax])
            bestz = self.xvect[imin+iFluxMin]
            merit = self.yvect[imin+iFluxMin]
        except Exception as e:
            print(e)
            bestz = -1
            merit = -1
            
        return bestz, merit
    
    def getZIndex(self, z):
        i=0
        for x in range(0,self.n-1):
            if z >= self.xvect[x]:
                i=x
            if z < self.xvect[x+1]:
                break
        return i
        
    def smoothListGaussian(self,list,degree=5):  
        window=degree*2-1  
        weight=np.array([1.0]*window)  
        weightGauss=[]  
        for i in range(window):  
            i=i-degree+1  
            frac=i/float(window)  
            gauss=1/(np.exp((4*(frac))**2))  
            weightGauss.append(gauss) 
        weight=np.array(weightGauss)*weight  
        smoothed=[0.0]*(len(list)-window)  
        for i in range(len(smoothed)):  
            smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight) 
        
        dec = len(list)-len(smoothed)
        decInit = int(dec/2.0)
        decEnd = dec-decInit
        #print dec
        smoothed2=[0.0]*(len(list)) 
        for i in range(0,decInit):
            smoothed2[i]=smoothed[0]
        for i in range(len(smoothed)):  
            smoothed2[i+decInit]=smoothed[i]
        for i in range(len(smoothed), len(list)-decInit):  
            smoothed2[i+decInit]=smoothed[len(smoothed)-1]
        return smoothed2 
     
    def getFluxExtrema(self, n):
        #x=np.array([6,3,5,2,1,4,9,7,8])
        #y=np.array([2,1,3,5,3,9,8,10,7])
        x = self.xvect
        y = self.smoothListGaussian(self.yvect,25)
        #y = self.yvect
        
        #print len(y)
        #print len(self.yvect)
        
        #pp.plot(self.xvect, self.yvect, 'k', label=self.name)
        #pp.plot(self.xvect, y, 'b', label=self.name)
        #pp.show()
        
        length = len(x)
        sortId=np.argsort(x)
        
        x=[x[s] for s in sortId]
        y=[y[s] for s in sortId]
        minm = np.array([])
        maxm = np.array([])
        i = 0
        while i < length-1:
            if i < length - 1:
                while i < length-1 and y[i+1] >= y[i]:
                    i+=1
        
                if i != 0 and i < length-1:
                    maxm = np.append(maxm,i)
        
                i+=1
        
            if i < length - 1:
                while i < length-1 and y[i+1] <= y[i]:
                    i+=1
        
                if i < length-1:
                    minm = np.append(minm,i)
                i+=1
 
 
        print("N min found : {0}".format(len(minm)))
        xmin = [self.xvect[int(s)] for s in minm]
        ymin = [self.yvect[int(s)] for s in minm]
        sortId=np.argsort(ymin)
        xmin2 = [xmin[s] for s in sortId]
        #ymin2 = [ymin[s] for s in sortId]
        for a in range(n):
            print("Min #{0} = {1}".format(a, xmin2[a]))
        return xmin2
       
    def getFluxExtrema2(self, n, enablePlot=False):
        #x=np.array([6,3,5,2,1,4,9,7,8])
        #y=np.array([2,1,3,5,3,9,8,10,7])
        xin = np.copy(self.xvect)
        yin = self.smoothListGaussian(self.yvect-self.getSplineContinuum(),25)
        #y = self.yvect
        
        #print len(y)
        #print len(self.yvect)

                
        length = len(xin)
        sortId=np.argsort(xin)
        
        x=[xin[s] for s in sortId]
        y=[yin[s] for s in sortId]
        minm = np.array([])
        maxm = np.array([])
        i = 0
        while i < length-1:
            if i < length - 1:
                while i < length-1 and y[i+1] >= y[i]:
                    i+=1
        
                if i != 0 and i < length-1:
                    maxm = np.append(maxm,i)
        
                i+=1
        
            if i < length - 1:
                while i < length-1 and y[i+1] <= y[i]:
                    i+=1
        
                if i < length-1:
                    minm = np.append(minm,i)
                i+=1
 
 
        print("getFluxExtrema2: N min found : {0}".format(len(minm)))
        xmin = [x[int(s)] for s in minm]
        ymin = [y[int(s)] for s in minm]
        sortId=np.argsort(ymin)
        xmin2 = [xmin[s] for s in sortId]
        #print sortId
        #print ymin
        #print xmin        
        #print xmin2
        #print("xmin2 = {}".format( xmin2)) 
        #ymin2 = [ymin[s] for s in sortId]
        for a in range(n):
            print("Min #{0} = {1}".format(a, xmin2[a]))
            
        if enablePlot:
            pp.plot(x, self.yvect, 'k', label=self.name)
            pp.plot(x, y, 'b', label=self.name)
            pp.show()
        
        return xmin2
       
    def getParabolicFitArea(self, z, dz=1e-4, constantLine=[],enablePlots = False):
        if enablePlots:
            print("getparabolicFit z: {}".format(z))
            
        if z<self.xvect[0] or z>self.xvect[len(self.xvect)-1]:
            return -1, z
           
        yprim = self.smoothListGaussian(self.yvect,10)
        if enablePlots:
            pp.plot(self.xvect, self.yvect, 'k', label=self.name)
            pp.plot(self.xvect, yprim, 'b', label=self.name)
            pp.show()
        
        # find indexes
        # Need to smooth before searching for the interval to fit ?
        for x in range(len(self.xvect)):
            if z<=self.xvect[x]:
                xz = x
                break

        for x in range(xz+1, len(self.xvect)):
            if yprim[x] < yprim[x-1]:
                xz = x
            else:
                break
        for x in range(xz-1, 0, -1):
            if yprim[x] < yprim[x+1]:
                xz = x
            else:
                break            
          
        xmax = len(self.xvect)-1
        for x in range(xz+1, len(self.xvect)):
            if yprim[x] < yprim[x-1]:
                xmax = x
                break
        xmin=0
        for x in range(xz-1, 0, -1):
            if yprim[x] < yprim[x+1]:
                xmin = x
                break
            
        xmin = int(xz - (xz-xmin)*1.0/4.0);
        xmax = int(xz + (xmax-xz)*1.0/4.0);
        
        if enablePlots:
            print("getparabolicFit interval: min = {0}, xz = {1}, xmax = {2}".format(xmin, xz, xmax))
        x1 = self.xvect[xmin:xmax]
        y1 = self.yvect[xmin:xmax]
        # Use polyfit.
        coefficients1 = np.polyfit(x1,y1,2)
        #print coefficients1
        if enablePlots:
            print("getparabolicFit coeff.: 0 = {0}, 1 = {1}, 2 = {2}".format(coefficients1[0], coefficients1[1], coefficients1[2]))
        #if coefficients1[0] < 0.8/dz:
        #    return -1, z
        
        polynomial = np.poly1d(coefficients1)
        # Feed data into pyplot.
        xpoints = np.linspace(self.xvect[xmin], self.xvect[xmax], xmax-xmin+1)
           
        parabole = self.xvect[xmin:xmax]
        droite = self.xvect[xmin:xmax]
        
        #if enablePlots:
        #        print("len(constantLine)={}, constantLine={}".format(len(constantLine), constantLine))
        if(len(constantLine)<2):
            a = (self.yvect[xmax]-self.yvect[xmin])/(self.xvect[xmax]-self.xvect[xmin])
            b = self.yvect[xmax] - a * self.xvect[xmax]
        else:
            a = constantLine[0]
            b = constantLine[1]
            if enablePlots:
                print("overriding base line for area detection: a={}, b={}".format(a, b))
                
        for x in range(len(x1)):
            #parabole[x] = np.polynomial.polynomial.polyval(x1[x],coefficients1)
            parabole[x] = coefficients1[0]*x1[x]*x1[x] + coefficients1[1]*x1[x] + coefficients1[2]
            droite[x] = a * x1[x] + b 
        
        if enablePlots:
            #pp.plot(x1,y1,'x',xpoints,polynomial(xpoints),'-')
            #print len(x1)
            #print len(parabole)
            pp.plot(x1,y1,'x', label = 'raw data')
            pp.plot(x1,parabole,'o', label = 'parabole fit')
            pp.plot(x1,droite,'-', label = 'continuum line')
            pp.legend()
            pp.xlabel('z')
            # Draw the plot to the screen
            pp.show()
        
        # estimate area between continuum line and parabole fit
        area = 0.0
        aminParabolicfit = 1e6
        for x in range(len(x1)):
            area += (droite[x] - parabole[x])/parabole[x]
            if aminParabolicfit > parabole[x]:
               aminParabolicfit = parabole[x]
               zminParabolicfit = x1[x]
               
        #override area by the analytic parabola area over a horiz. line: 2/3 * b * h
        print("coefficients1[2]*dz = {}, aminParabolicfit={}".format(coefficients1[2]*dz, aminParabolicfit ))
        H = -(aminParabolicfit-b)       
        if H<0:
            H=0
        B = 2 * math.sqrt(abs(H/(coefficients1[0]*dz)))
        
        area = (2.0/3.0) * B * H                       
               
        #if enablePlots:
        print("getParabolicFit equ. area: area= {0}".format(area))
            
        
        
        
        return area, zminParabolicfit

    def get_StrongELSnrPrior(self):
        coeff = 10.0
        priors = [a*coeff for a in self.amazed_StrongELSNR]
        return priors
        
    def getMeritExtremum(self, idxExtremum):
        z = self.amazed_extrema[idxExtremum]
        print("searching merit for extrema_idx = {}, z={}".format(idxExtremum, z))
        ind = np.argmin(np.abs(np.array(self.xvect)-z))
        merit = self.yvect[ind]
        
        print("Found merit for chi2_idx = {}, merit={}".format(ind, merit))
        return merit
        
            

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./chisquare.py -s """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-i", u"--input", help="path to the chi2 curve to be plotted",  dest="chi2Path", default="")
    parser.add_option(u"-o", u"--input2", help="path to the other chi2 curve to be plotted",  dest="chi2Path2", default="")

    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        print('using full path: {0}'.format(options.chi2Path))
        s = ResultChisquare(options.chi2Path)
        print(s) 
        
        if options.chi2Path2 == "":    
            s.plot(showContinuumEstimate=False, showExtrema=True, showAmbiguities=False)
        else:
            s2 = ResultChisquare(options.chi2Path2) 
            s.plotCompare(s2) 
    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()
        
 
if __name__ == '__main__':
    print("Chisquare")
    if 1:
        Main( sys.argv )
        
      
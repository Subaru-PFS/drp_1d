# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 14:30:52 2015

@author: aschmitt
"""

import sys
import os
import time
import optparse
import matplotlib.pyplot as plt
import numpy as np
import math

import resparser as rp
import templatechi as tchi
import resultstat as rstat
import chisquare as chisq                
                
class BestRedshift(object):
    def __init__(self, dir, spcname, zlist, method = "spctplchi2sarea", chi2type="linemodel"):
        self.logTagStr = "BestRedshift"
        self.dir = dir
        self.spcname = spcname
        self.zlist = zlist
        #self.type = "spctplchi2sweighted"
        #self.type = "spctplchi2sarea"
        self.type = method
        
        #self.zlist = [1.3455, 2.1634] #atm clean 107, chi2nocontinuum gives correct z, the tplprior doesn't help
        #self.zlist = [1.1777, 0.4178] #spc_example_sc_020089640_F02P019_vmM1_red_93_1_atm_clean, chi2nocontinuum could help,Z emission vs Z galaxy, tplprior doesn't help 
        #self.zlist = [ 1.134, 0.8876] #sc_020106186_F02P018_vmM1_red_18_1_atm_clean, tplprior doesn't help
        #self.zlist = [ 1.3698, 1.125] #sc_020115875_F02P019_vmM1_red_33_1_atm_clean, 2nd minimum in the chi2, tplprior doesn't help
        #self.zlist = [ 1.1654, 1.5851] #sc_020132417_F02P017_vmM1_red_93_1_atm_clean, OK with tplprior
        #self.zlist = [ 1.3202, 0.2928] #sc_020141834_F02P030_vmM1_red_130_1_atm_clean, 2nd minimum in the chi2nocontinuum, not all templates at z=1.32, tplprior doesn't help
        #self.zlist = [ 1.0322, 0.8892] #sc_020142252_F02P031_vmM1_red_98_1_atm_clean
        #self.zlist = [ 1.3052, 0.636] #sc_020144757_F02P029_vmM1_red_37_1_atm_clean, not all templates at z=1.30, tplprior doesn't help much
        
        # weighting = coeff = 31/97 corrected
        # weighting = coeff, when tplchival[i] forced to 0 = 18/97 corrected
        # weighting = coeff, when limited to exact same support = 34/97 corrected
        
        # uing an automatic z extrema search with smoothing, 3 extrema : 5/97 corrected
        
        print("#self.zlist = [ {0}, {1}] #{2},".format(self.zlist[0], self.zlist[1], self.spcname))  
        
        self.spctplchi2s = None
        self.spctplchi2sweighted = None
        self.tplchi2s = None
        self.spctplchi2sarea = None
        self.spctplchi2szcenter = None
        
        self.matrix = None
        
        self.tplPlotk = 2
        self.tplPlotWeight = 1.0
        
        
        if self.type=="spctplchi2sweighted":
            self.loadTplChi2Weighting(chi2type)
        elif self.type=="spctplchi2sarea":
            self.loadArea(chi2type)
        
        
    def loadTplChi2Weighting(self, chi2type="linemodel"):
        s = rp.ResParser(self.dir)
        print(s) 
        
        # load spc tpl chi2s
        for x in range(0,len(self.zlist)):
            chi2list = s.getChi2List(self.spcname, self.zlist[x], chi2type=chi2type)
            print('{0} spctplchi load : {1}'.format(x, chi2list))
            if x==0:
                self.spctplchi2s = np.zeros((len(self.zlist), len(chi2list)))
                self.spctplchi2sweighted = np.zeros((len(self.zlist), len(chi2list)))
            for k in range(len(chi2list)):
                self.spctplchi2s[x,k] = chi2list[k]
                self.spctplchi2sweighted[x,k] = chi2list[k]
        print('{0} spctplchis loaded.'.format(len(self.zlist)))

        # load tpl chi2s
        tplDir = s.getConfigVal('templatedir')
        tplchi = tchi.TemplateChi(tplDir)
        tplchi.importTxt()
        self.tplchi2s = tplchi.tplchi2s
        print('{0} tplchi loaded : {1}'.format(1, self.tplchi2s.shape))
        
    def loadArea(self, chi2type="linemodel"):
        s = rp.ResParser(self.dir)
        print(s) 
        nmodels = 1
        tplpaths = [""]
        if  not chi2type=="linemodel":
            tplpaths = s.getTplFullPathList()
            nmodels = len(tplpaths)

        # load 
        for k in range(nmodels):
            filepath = s.getChi2FullPath(self.spcname, os.path.basename(tplpaths[k]), chi2type=chi2type)
            print("loading area for tpl: {0}".format((k)))
            chi2 = chisq.ResultChisquare(filepath) 
            for x in range(0,len(self.zlist)):
                if  chi2type=="linemodel":
                    hLine = [0.0,chi2.getFluxMax()]
                else:
                    hLine = []
                area, zcenter = chi2.getParabolicFitArea(self.zlist[x], dz=chi2.xvect[1]-chi2.xvect[0], constantLine=hLine, enablePlots=False)
            
                if x==0 and k==0:
                    self.spctplchi2sarea = np.zeros((len(self.zlist), len(tplpaths)))
                    self.spctplchi2szcenter = np.zeros((len(self.zlist), len(tplpaths)))
                self.spctplchi2sarea[x,k] = area
                self.spctplchi2szcenter[x,k] = zcenter
                #print("saved {0} in location {1},{2}".format(area, x,k))
        print('{0} spctplchisareas loaded.'.format(len(self.zlist)))


    def applyWeighting(self):
        print("Applying weight...")
        for y in range(self.spctplchi2s.shape[1]):
            #find the commom support
            indValid = np.ones(self.spctplchi2s.shape[1])
            for x in range(self.spctplchi2s.shape[0]):
                for itpl in range(self.spctplchi2s.shape[1]):
                    if self.spctplchi2s[x,itpl]>1e6:
                        indValid[itpl] = 0
            for x in range(self.spctplchi2s.shape[0]):
                w,a = self.weight(self.spctplchi2s[x,:], self.tplchi2s[y,:], y, indValid)
                #a = 10.0
                #w = self.weightgetVal(self.spctplchi2s[x,:], self.tplchi2s[y,:], a)
                
                print("weight min found {0}, for z = {1} : a = {2}, w = {3}".format(x,self.zlist[x],a,w))
                coeff = w#math.log(w)
                if x==0 and y==self.tplPlotk:
                    self.tplPlotWeight = a*coeff
                self.spctplchi2sweighted[x,y] = self.spctplchi2s[x,y]*coeff
        
    def weight(self, spc, s2, iy, indValid):
#        n=len(spc)
#        sumspc = 0.0
#        for x in range(n):            
#            if(spc[x]<1e6):
#                sumspc+=spc[x]
#        for x in range(n):            
#            if(spc[x]<1e6):
#                spc[x] -= sumspc/n 
#        n=len(s2)
#        sumspc = 0.0
#        for x in range(n):            
#            if(s2[x]<1e6):
#                sumspc+=s2[x]
#        for x in range(n):            
#            if(spc[x]<1e6):
#                s2[x] -= sumspc/n 
#              
        # translate to make 0 on the cirrent template row
        sc1 = []
        sc2 = []
        n=len(spc)
        #if(abs(spc[iy])<1e6):
        #    print("{0} spc = {1}".format(iy, spc[iy]))
        for x in range(n): 
            if 0:
                sc1.append(spc[x] - spc[iy])
            else:
                sc1.append(spc[x])   
            #if(abs(spc[iy])<1e6):  
            #    sc1.append(spc[x] - spc[iy])
            #else:
            #    sc1.append(spc[x])
        for x in range(n):  
            if 0:
                sc2.append(s2[x] - s2[iy])
            else:
                sc2.append(s2[x])       
            #if(abs(s2[iy])<1e6):  
            #    sc2.append(s2[x] - s2[iy])
            #else:
            #    sc2.append(s2[x])
            
        a = self.weightgetA(sc1, sc2, indValid)
        #print("a = {0}".format(a))
        na = 1000
        amin = a/100
        amax = a*100
        astep = (amax-amin)/na
        aylist = []
        axlist = []
        for x in range(na):
            A = amin + astep*x
            aylist.append(self.weightgetVal(sc1, sc2, A, indValid))
            axlist.append(A)
            
        if 0:
            plt.plot(axlist, aylist)
            plt.show()
        
        sortId=np.argsort(aylist)
        ay=[aylist[s] for s in sortId]
        ax=[axlist[s] for s in sortId]
        return ay[0],ax[0]
        
    def weightgetA(self, spc, s2, indValid):
        #estimate A
        n=len(spc)
        sumsq1 = 0.0
        sumsq2 = 0.0
        for x in range(n):
            if(spc[x]<1e6 and s2[x]<1e6) and indValid[x]==1:
                sumsq1 += spc[x]
                sumsq2 += s2[x]  
        A = abs(sumsq1/sumsq2)
        #A = 1.0
        return A
        
    def weightgetVal(self, spc, s2, A, indValid):
        sum2=0.0
        n = len(spc)
        nsum = 0.0
        for x in range(n):
            if(spc[x]<1e6 and s2[x]<1e6) and indValid[x]==1:
                a = spc[x]-A*s2[x]
                asq = a*a
                sum2 += asq  
                nsum += 1
        #print("nsum = {0}".format(nsum))
        return sum2/nsum

    def plot(self, titletag=""):
        if self.type=="spctplchi2sweighted":
            self.matrix = self.spctplchi2sweighted 
            basePlotName = "spctplchi2"
            yName = "chi2"
        elif self.type=="spctplchi2sarea":
            self.matrix = self.spctplchi2sarea
            basePlotName = "spctplchi2area"
            yName = "chi2 area"
                
                
         #find limits
        cmin = +1e6
        cmax = -1e6
        thres = 1e6
        legendz = []
        for x in range(self.matrix.shape[0]):
            legendz.append(self.zlist[x])
            for y in range(self.matrix.shape[1]):
                if cmin >  self.matrix[x,y] and self.matrix[x,y]>-thres:
                    cmin = self.matrix[x,y]
                if cmax <  self.matrix[x,y] and self.matrix[x,y]<thres:
                    cmax = self.matrix[x,y]
                  
        plt.figure("GetBestRedshift_{0}".format(titletag))
        plt.plot(np.transpose(self.matrix))
        plt.xlabel('template index')
        plt.ylabel('{0}'.format(yName))
        name1 = "{0} for {1} \n{2}".format(basePlotName, self.spcname, titletag)
        plt.title(name1)
        
        if self.type=="spctplchi2sweighted":
            tplPloty = []
            for x in range(len(self.tplchi2s[self.tplPlotk,:])):
                tplPloty.append(self.tplPlotWeight*self.tplchi2s[self.tplPlotk,x])
            plt.plot(tplPloty)
            legendz.append("tpl")
        plt.legend(legendz)
        plt.grid()
        plt.ylim([cmin,cmax])
        plt.show()
        
    def matplot(self, titletag=""):
        if self.type=="spctplchi2sweighted":
            self.matrix = self.spctplchi2sweighted 
        elif self.type=="spctplchi2sarea":
            self.matrix = self.spctplchi2sarea
            
         #find limits
        cmin = +1e6
        cmax = -1e6
        thres = 1e6
        legendz = []
        for x in range(self.matrix.shape[0]):
            legendz.append(self.zlist[x])
            for y in range(self.matrix.shape[1]):
                if cmin >  self.matrix[x,y] and self.matrix[x,y]>-thres:
                    cmin = self.matrix[x,y]
                if cmax <  self.matrix[x,y] and self.matrix[x,y]<thres:
                    cmax = self.matrix[x,y]
                  
        #plt.figure("matplot".format(titletag))
        #plt.plot(np.transpose(self.matrix))
        plt.matshow(np.transpose(self.matrix), interpolation='nearest',aspect='auto')
        
        plt.xlabel('z')
        plt.ylabel('template')
        name1 = "spctplchi2 for {0} \n{1}".format(self.spcname, titletag)
        plt.title(name1)
        

        #plt.legend(legendz)
        #plt.grid()
        plt.colorbar()
        plt.clim(cmin,cmax)
        plt.clim(cmin,cmin*2.0)
        plt.show()
        
    def getBestRedshift(self):
        if self.type=="spctplchi2sweighted":
            self.matrix = self.spctplchi2sweighted 
        elif self.type=="spctplchi2sarea":
            self.matrix = self.spctplchi2sarea
            
         #find limits
        cmin = +1e6
        tmin = [0,0]
        cmax = -1e6
        tmax = [0,0]
        thres = 1e6
        legendz = []
        for x in range(self.matrix.shape[0]):
            legendz.append(self.zlist[x])
            for y in range(self.matrix.shape[1]):
                if cmin >  self.matrix[x,y] and self.matrix[x,y]>-thres:# and y!=16:
                    cmin = self.matrix[x,y]
                    tmin = [x,y]
                if cmax <  self.matrix[x,y] and self.matrix[x,y]<thres:
                    cmax = self.matrix[x,y]
                    tmax = [x,y]
        print('Max L : {}, best redshift = {} found for x(=i_z) = {} and y(=i_tpl) = {}'.format(cmax, self.zlist[tmax[0]], tmax[0], tmax[1]))
        print('Min L : {}, best redshift = {} found for x(=i_z) = {} and y(=i_tpl) = {}'.format(cmin, self.zlist[tmin[0]], tmin[0], tmin[1])) 
        
        if self.type=="spctplchi2sweighted":
            return tmin[0] 
        elif self.type=="spctplchi2sarea":
            print("best redshift = {0}, corrected by parabolic fit = {1}".format(self.zlist[tmax[0]], self.spctplchi2szcenter[tmax[0], tmax[1]]))
            return tmax[0]
        
  
#def getZCandidatesFromDiff(resList, indice=0):
#    return [resList.list[indice].zref, resList.list[indice].zcalc]
#
#def getZCandidatesFromChi2Extrema(resDir, resList, indice=0):
#    redshiftslist = []
#    spcName = resList.list[indice].name
#    s = rp.ResParser(resDir)
#    tplpaths = s.getTplFullPathList()
#    for a in range(len(tplpaths)):
#        filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), "nocontinuum")
#        chi2 = chisq.ResultChisquare(filepath)
#        nextrema = 3
#        redshiftstpl = chi2.getFluxExtrema(nextrema)
#        for b in range(nextrema):
#            redshiftslist.append(redshiftstpl[b])
#    redshiftsUnsorted = list(set(redshiftslist))
#    sortId=np.argsort(redshiftsUnsorted)
#    redshifts = [redshiftsUnsorted[b] for b in sortId]
#    return redshifts
#
#def getZCandidatesFromAmazedChi2Extrema(resDir, resList, indice=0):
#    redshiftslist = []
#    spcName = resList.list[indice].name
#    s = rp.ResParser(resDir)
#    tplpaths = s.getTplFullPathList()
#    for a in range(len(tplpaths)):
#        #filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), "nocontinuum")
#        filepath = s.getChi2FullPath(spcName, os.path.basename(tplpaths[a]), "raw")
#        chi2 = chisq.ResultChisquare(filepath)
#        nextrema = 5
#        redshiftstpl = chi2.amazed_extrema
#        for b in range(nextrema):
#            redshiftslist.append(redshiftstpl[b])
#    redshiftsUnsorted = list(set(redshiftslist))
#    sortId=np.argsort(redshiftsUnsorted)
#    redshifts = [redshiftsUnsorted[b] for b in sortId]
#    return redshifts
#    
        
def do(resDir, spcName, opt=0):
    print("\n\n")    
    print('using amazed results full path: {0}'.format(resDir))
    print('processing spc : {0}'.format(spcName))
    
    chi2Type = "raw"
    
    print("\n")    
    resList = rstat.ResultList(resDir, diffthreshold=0.01, opt='brief')        
    print("Results, N found failures = {0}".format(resList.n))
    if spcName == "":
        indice = 3
        if 1:
            redshifts = resList.getZCandidatesFromDiff(indice)
            #redshifts = resList.getZCandidatesFromChi2Extrema(indice)
            #redshifts = getZCandidatesFromChi2Extrema(resDir, resList, indice)
            #redshiftsUnsorted = [0.398900, 2.090900, 0.422300, 1.169900, 0.404300, 1.585100, 1.165400, 1.165000, 0.983500, 1.060900, 1.165700, 1.165000, 1.163900, 0.421500, 0.406800, 0.421200, 0.409200]    
            #redshiftsUnsorted = [0.398900,1.585100, 1.165400,  1.060900]   
            #redshiftsUnsorted = [1.363, 1.3808 , 2.9231] 
            #sortId=np.argsort(redshiftsUnsorted)
            #redshifts = [redshiftsUnsorted[b] for b in sortId]
        else:
            redshifts, m = resList.getZCandidatesFromAmazedChi2Extrema(indice = indice, chi2Type = chi2Type, nextrema = 10)
            #redshifts.append(1.3455)
            #sortId=np.argsort(redshifts)
            #redshifts = [redshifts[b] for b in sortId]
        

        spcName = resList.list[indice].name
    else:
        for x in range(resList.n):
            if spcName == resList.list[x].name:
                indice = x
                break
        #redshifts = resList.getZCandidatesFromDiff(indice)
        redshifts, m = resList.getZCandidatesFromAmazedChi2Extrema(indice=indice, chi2Type = "linemodel", nextrema = 5)
        #redshiftsUnsorted = [1.363, 1.3808 , 2.9231] 
        #sortId=np.argsort(redshiftsUnsorted)
        #redshifts = [redshiftsUnsorted[b] for b in sortId]
    print("Redshift candidates: {0}".format(redshifts))

    print("Processing spectrum #{0} : {1}".format(indice, spcName))
    #method = "spctplchi2sarea"
    method = "spctplchi2sweighted"
    bestz = BestRedshift(resDir, spcName, redshifts, method, chi2Type)  
    
    if method == "spctplchi2sarea" :  
        z0 = redshifts[bestz.getBestRedshift()]
        print("\n")
        print("{0} : zref = {1}, zarea = {2}".format(spcName, resList.list[indice].zref, z0))
        print("\n")        
        bestz.plot("AREA")
    else:   
        z0 = redshifts[bestz.getBestRedshift()]
        bestz.plot("RAW")
        #bestz.applyWeighting()
        #z1 = redshifts[bestz.getBestRedshift()]
        print("\n")
        #print("{0} : zref = {1}, zweighted = {2}".format(spcName, resList.list[indice].zref, z1))
        print("{0} : zref = {1}, zraw = {2}".format(spcName, resList.list[indice].zref, z0))
        print("\n")
        #bestz.plot("WEIGHTED")
    
 
def doAll(resDir, opt=0):
    print("\n\n")    
    print('using amazed results full path: {0}'.format(resDir))
    print('processing all spc ')
    
    print("\n")    
    resList = rstat.ResultList(resDir, diffthreshold=0.01, opt='brief')        
    print("Results, N found failures = {0}".format(resList.n))
    correctedN = 0
    for x in range(0,resList.n):
        spcName = resList.list[x].name
        indice = x
        if 0: # redshifts from stat file
            redshifts = resList.getZCandidatesFromDiff(indice)
        else: #redshift from chi2nocontinuum extrema
            #redshifts = resList.getZCandidatesFromChi2Extrema(indice)
            redshifts = resList.getZCandidatesFromAmazedChi2Extrema(indice)
        print("Redshift candidates: {0}".format(redshifts))
        
        print("\n\n")  
        print("Processing spectrum #{0} : {1}".format(indice, spcName))    
        bestz = BestRedshift(resDir, spcName, redshifts)    
        z0 = redshifts[bestz.getBestRedshift()]
        #bestz.plot("RAW")
        bestz.applyWeighting()
        z1 = redshifts[bestz.getBestRedshift()]
        #bestz.plot("WEIGHTED")  
        if abs(z1-resList.list[indice].zref)<0.01:
            correctedN+=1
        print("{0} : zref = {1}, zweighted = {2}".format(spcName, resList.list[indice].zref, z1))
        print("{0} : zref = {1}, zraw = {2}".format(spcName, resList.list[indice].zref, z0))        
        print("Failures corrected = {0}/{1}".format(correctedN, x+1))
        time.sleep(1)
        #raw_input("PRESS ENTER TO CONTINUE.")
          
    print("GetBestRedshift WEIGHTING STATS: Failures corrected = {0}/{1}".format(correctedN, resList.n))

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./getbestredshift.py """
    parser = optparse.OptionParser(usage=usage)
    #parser.add_option(u"-x", u"--xaxis", help="force the use of indexes for the xaxis (--xaxis=""index""), or enable the use of calibrated wavelength (--xaxis=""angstrom"")",  dest="xAxisType", default="angstrom")
    parser.add_option(u"-d", u"--dir", help="path to the amazed results directory (/output/)",  dest="resDir", default="./output")
    parser.add_option(u"-s", u"--spc", help="name of the spectrum to be plotted",  dest="spcName", default="")
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        do(options.resDir, options.spcName)
        #doAll(options.resDir)
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
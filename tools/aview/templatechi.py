# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:21:16 2015

@author: aschmitt
"""
import sys
import os
import optparse

import matplotlib as mpl
mpl.use('Qt5Agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import scipy.stats as sps

import resparser as rp
import spectrum as sp

class TemplateChi(object):
    def __init__(self, dirpath):
        self.logTagStr = "TemplateChi"
        self.path = dirpath
        self.ntpl = -1.0;
        self.tplPathList = []
        self.tplchi2s = None
        # template categories
        self.config_tplCategories = ['emission', 'galaxy', 'qso', 'star']
        self.lambdarange = [3000, 9600]
        self.load()
                
    def load(self):
        self.tplPathList = self.getTplFullPathList()
        self.ntpl = len(self.tplPathList)
        #self.ntpl = 3
        for x in range(0,self.ntpl):
            print("{0} : {1}".format(x, self.tplPathList[x]))

    def exportTxt(self):
        path = os.path.join(self.path, "tplchi2")
        np.savetxt(path, self.tplchi2s)

    def importTxt(self):
        pathroot = '/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2'
        #path = os.path.join(pathroot, "tplchi2_20150810csym")
        path = os.path.join(pathroot, "tplchi2_20150812e")
        self.tplchi2s = np.loadtxt(path, self.tplchi2s)

    def getTplFullPathList(self):
        """
        """
        tplPathList = []
        tplRootPath = self.path
        for cat in self.config_tplCategories:
            tplPath = os.path.join(tplRootPath, cat)
            if os.path.exists(tplPath):
                for file in sorted(os.listdir(tplPath)):
                    tplPathList.append(os.path.join(tplPath, file)) 
        return tplPathList
        
    def computeChi2Ref(self):
        """
        """
        self.tplchi2s = np.zeros((self.ntpl, self.ntpl))
        for x in range(self.ntpl):
            t1 = sp.Spectrum(self.tplPathList[x], 'template', True)
            #for y in range(0, x):
            #    chi2 = self.tplchi2s[y,x]        
            #    print("chi2_ {0},{1} = {2} : SYM".format(x, y, chi2))  
            #    self.tplchi2s[x,y]=chi2
            for y in range(0, self.ntpl):
                t2 = sp.Spectrum(self.tplPathList[y], 'template', True)
                #chi2 = sps.chisquare(t.yvect, t.yvect)
                t2b,n = self.rebin(t1, t2)
                #plt.plot(t2.xvect, t2.yvect, 'ko') 
                #plt.plot(t2b.xvect, t2b.yvect, 'ro') 
                #plt.show()        
                chi2 = self.chi2(t1, t2b, n)        
                print("chi2_ {0},{1} = {2}".format(x, y, chi2))  
                self.tplchi2s[x,y]=chi2
        #plt.plot(chi2) 
        #plt.show()
        self.exportTxt()

    def rebin(self, s1, s2):
        s2b = s1.copy()
        n=0
        for x in range(len(s1.xvect)):
            if s1.xvect[x] >= s2.xvect[0] and s1.xvect[x] <= s2.xvect[len(s2.xvect)-1] and s1.xvect[x] > self.lambdarange[0] and s1.xvect[x] < self.lambdarange[1]: 
                s2b.yvect[x] = np.interp(s1.xvect[x], s2.xvect, s2.yvect)
                n +=1
            else:
                s2b.yvect[x] = s1.yvect[x]
        return s2b, n

    def chi2(self, s1, s2, ntrue):
        #estimate A
        sq1 = [a*a for a in s1.yvect]
        sumsq1 = sum(sq1)
        sq2 = [a*a for a in s2.yvect]
        sumsq2 = sum(sq2)
        A = sumsq1/sumsq2
        sum2=0.0
        n = len(s1.xvect)
        for x in range(n):
            a = s1.yvect[x]-A*s2.yvect[x]
            asq = a*a
            sum2 += asq    
        #return A
        return sum2/ntrue
        
    def plot(self, k):
        print(self.tplchi2s[:,k])
        
        ymin = min(self.tplchi2s[:,k])
        ymax = max(self.tplchi2s[:,k])
        if ymax>1.0:
            ymax = 1.0
        
        #self.importTxt()
        plt.plot(self.tplchi2s[:,k])
        plt.xlabel('index')
        plt.ylabel('chi2')
        name1 = "tplchi2 for {0} \nindex = {1}".format(os.path.basename(self.tplPathList[k]), k)
        plt.title(name1)
        plt.ylim([ymin,ymax])
        plt.show()
        
    def matplot(self):
        #find limits
        cmin = +1e6;
        cmax = -1e6;
        thres = 1e6
        for x in range(self.ntpl):
            for y in range(self.ntpl):
                if cmin >  self.tplchi2s[x,y] and self.tplchi2s[x,y]>-thres:
                    cmin = self.tplchi2s[x,y]
                if cmax <  self.tplchi2s[x,y] and self.tplchi2s[x,y]<thres:
                    cmax = self.tplchi2s[x,y]
                
        print("cmin={0}, cmax={1}".format(cmin, cmax))
        #fig = plt.figure()
        #ax = fig.add_subplot(1,1,1)
        #ax.set_aspect('equal')
        plt.matshow(self.tplchi2s, interpolation='nearest', cmap=plt.cm.ocean)
        plt.colorbar()
        plt.clim(cmin,cmax)
        #plt.clim(cmin,2)
        plt.show()
        
        
def do(tplDir, opt=0):
    print('using full path: {0}'.format(tplDir))
    tplchi = TemplateChi(tplDir)
    #tplchi.computeChi2Ref()
    
    tplchi.importTxt()
    tplchi.plot(0)
    tplchi.matplot()

    

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./templatechi.py """
    parser = optparse.OptionParser(usage=usage)
    #parser.add_option(u"-x", u"--xaxis", help="force the use of indexes for the xaxis (--xaxis=""index""), or enable the use of calibrated wavelength (--xaxis=""angstrom"")",  dest="xAxisType", default="angstrom")
    parser.add_option(u"-d", u"--dir", help="path to the template catalog directory (/Templates/)",  dest="tplDir", default="./Templates")
    #parser.add_option(u"-f", u"--diff", help="path to the diff file",  dest="diffFile", default="diff.txt")
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
	do(options.tplDir)
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
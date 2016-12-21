# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 20:12:38 2016

@author: aschmitt
"""

import os
import sys

import argparse

import reference

class Spectrumlist(object):
    def __init__(self, path):
        self.logTagStr = "SpcList"
        
        self.path = path
        self.name = os.path.basename(path)

        self.fvect = []
        self.errfvect = []
            
        self.load()
 

    def load(self):
        """
        load the spectrumlist file data
        """ 
        print("load".format())

        f = open(self.path)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                data = lineStr.split("\t")
                #if len(data)<2:
                #    data = lineStr.split(" ")
                #data = [r for r in data if r != '']
                self.fvect.append(data[0])
                if len(data)>=2:
                    self.errfvect.append(data[1])
        print("Spctrumlist loaded n F = {}".format(len(self.fvect)))
        print("Spctrumlist loaded n ErrF = {}".format(len(self.errfvect)))
        
              
        
    def filterSpectrumList(self, refFilePath, magmin=-1, magmax=50, zmin=-1, zmax=50, sfrmin=-1, sfrmax=1e6, outputPath=""):
        ref = reference.Reference(refFilePath)
        idList = self.fvect
        indexes = ref.filterIdList(idList, magmin=magmin, magmax=magmax, zmin=zmin, zmax=zmax, sfrmin=sfrmin, sfrmax=sfrmax)
        print("spllist filtering: indexes found n = {}".format(len(indexes)))

        if outputPath=="":
            dirPath = os.path.split(self.path)[0]
            nameNoExt = os.path.splitext(self.name)[0]
            outputFileFullPath = os.path.join(dirPath, "{}_z{}-{}_mag{}-{}_sfr{}-{}.spectrumlist".format(nameNoExt, zmin, zmax, magmin, magmax, sfrmin, sfrmax))
            
            f = open(outputFileFullPath, 'w')
            for k, idThis in enumerate(self.fvect):
                if k in indexes:
                    if len(self.errfvect) == len(self.fvect):
                        f.write("{}\t{}\n".format(self.fvect[k], self.errfvect[k]))
                    else:                        
                        f.write("{}\n".format(self.fvect[k]))
                        
    def splitIntoSubsets(self, subsetCount, outputDirPath):
        """
        Split the spectrumlist into N subsets with subsetCount items in each subsets
        """
        verbose = 0
        print("INFO: writing subsets with count: {}".format(subsetCount))
        print("INFO: writing subsets files in directory: {}".format(outputDirPath))
        iCountSubset = 0
        iCountInSubset = 0
        for k, f in enumerate(self.fvect): 
            if verbose:                        
                print("INFO: k={}".format(k))
                print("INFO: iCountInSubset={}, subsetCount={}".format(iCountInSubset, subsetCount))
                print("type iCountInSubset={}".format(type(iCountInSubset)))
                print("type subsetCount={}".format(type(subsetCount)))
            if iCountInSubset==subsetCount:
                if verbose:
                    print("INFO: switching subset file, iCountSubset={}".format(iCountSubset))
                fsub.close()
                iCountInSubset = 0
                iCountSubset+=1
                if verbose:
                    print("INFO: switching subset file, after inc. iCountSubset={}".format(iCountSubset))
                
            if iCountInSubset==0:
                name_noext = os.path.splitext(self.name)[0]
                subsetFilename = "{}_{}.spectrumlist".format(name_noext, iCountSubset)
                subsetFilePath = os.path.join(outputDirPath, subsetFilename)
                fsub = open(subsetFilePath, 'w')
                print("INFO: opening subset file: {}".format(subsetFilename))
                
            if len(self.errfvect) == len(self.fvect):
                fsub.write("{}\t{}\n".format(self.fvect[k], self.errfvect[k]))
            else:                        
                fsub.write("{}\n".format(self.fvect[k]))
            
            if verbose:
                print("INFO: inc. iCountInSubset={}".format(iCountInSubset))
            iCountInSubset+=1
            if verbose:
                print("INFO: after inc. iCountInSubset={}".format(iCountInSubset))
        
        try:
            fsub.close()
        except:
            pass
        
                        
        
        
def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--listPath", dest="listPath", default="",
                    help="path to the list file to be loaded") 
                    
    #options to filter the spectrumlist
    parser.add_argument("-f", "--filter", dest="filter", action='store_true',
                    help="enable spectrumlist filtering")
    parser.add_argument("-r", "--refPath", dest="refPath", default="",
                    help="path to the reference file to be loaded")    
                    
    parser.add_argument("-m", "--magRange", dest="magRange", default="-1.0 50.0",
                    help="magnitude range filter for the histograms")
    parser.add_argument("-s", "--sfrRange", dest="sfrRange", default="-1.0 10000.0",
                    help="sfr range filter for the histograms")
    parser.add_argument("-z", "--zRange", dest="zRange", default="-1.0 20.0",
                    help="redshift range filter for the histograms")   
                    
    #option to split the spectrumlist
    
    parser.add_argument("-d", "--divide", dest="divide", action='store_true',
                    help="enable spectrumlist dividing/splitting into subsets")
    parser.add_argument("-n", "--dividecount", dest="dividecount", default=100,
                    help="dividing/splitting count") 
          
    options = parser.parse_args()
    #print(options)

    if os.path.exists(options.listPath) :
        rpath = os.path.abspath(options.listPath)
        print('using full path: {}'.format(rpath))
        spclist = Spectrumlist(rpath)
 
        if options.filter:
            if not os.path.exists(options.refPath):
                print("ERROR: ref path does not exit, aborting...".format(options.refPath))
            else:
                print('filtering using ref. path: {}'.format(options.refPath))
            zRange = [-1.0, 20.0]
            zRange[0] = float(options.zRange.split(" ")[0])
            zRange[1] = float(options.zRange.split(" ")[1])
            magRange = [-1.0, 50.0]
            magRange[0] = float(options.magRange.split(" ")[0])
            magRange[1] = float(options.magRange.split(" ")[1])
            sfrRange = [-1.0, 10000.0]
            sfrRange[0] = float(options.sfrRange.split(" ")[0])
            sfrRange[1] = float(options.sfrRange.split(" ")[1]) 
            
            spclist.filterSpectrumList(options.refPath, magmin=magRange[0], magmax=magRange[1], zmin=zRange[0], zmax=zRange[1], sfrmin=sfrRange[0], sfrmax=sfrRange[1])
        elif options.divide :
            dirpath = os.path.split(options.listPath)[0]
            rpath = os.path.abspath(dirpath)
            outputPath = os.path.join(rpath, "spectrumlist_subs_{}".format(options.dividecount))
            if not os.path.exists(outputPath):
                os.mkdir(outputPath)
            print('splitting using full path: {}'.format(outputPath))
            spclist.splitIntoSubsets(int(options.dividecount), outputPath)
    else :
        print("Error: invalid arguments")
        exit()
    
    
    
def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print("Spectrumlist")
    Main( sys.argv )
    

# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 13:32:17 2017

@author: aschmitt
"""

### python 3 compatibility imports ###
#from __future__ import division, print_function
#from future_builtins import *
#from builtins import range
#from future.builtins.disabled import apply, cmp, file, raw_input, xrange
### ##### # ############# ####### ###

import os
import sys
import time
import argparse
import shutil



class processRecombine(object):
    def __init__(self, path, outputPath):
        self.logTagStr = "processRecombine"
        self.outputPath = outputPath
        self.flist = self.listInfoFiles(path)
        self.recombineAll()
  
    def listInfoFiles(self, directory, exclstring="xrfctgrefardtasrdtrdsatrdtartcomplicatedstring"):
        """
        Scans input directory for the .info files listing the subset directories
        """
        fileList = []
    	
        for i,file in enumerate(sorted(os.listdir(directory))):
            if file.endswith(".info") and exclstring not in file:  
                a1 = str(file)
                #print("file found ={}".format(a1))
                fullPath = os.path.join(directory, a1)
                fileList.append(fullPath)
        return fileList
        
    def recombine(self, infoFilePath):
        """
        recombine from the subset described by infoFilePath
        """
        subpathsList = []
        #load info file
        f = open(infoFilePath, 'r')
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                subpathsList.append(lineStr)
            
        #create this config merged output path
        name = os.path.split(infoFilePath)[1]
        name_noext = os.path.splitext(name)[0]
        outputPath = os.path.join(self.outputPath, name_noext)
        if not os.path.exists(outputPath):
            os.mkdir(outputPath)
            
        #merge redshift.csv files
        self.mergeCsvFiles(subpathsList, outputPath, fileName="redshift.csv")
        
        #merge log files
        self.mergeCsvFiles(subpathsList, outputPath, fileName="log.txt")

        try:        
            #merge input spectrumlist files
            self.mergeCsvFiles(subpathsList, outputPath, fileName="input.spectrumlist")         
            
            #copy json file from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "parameters.json")
            dest_path = os.path.join(outputPath, "parameters.json")
            shutil.copy(source_path, dest_path)
            
            #copy linecatalog file from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "linecatalog.txt")
            dest_path = os.path.join(outputPath, "linecatalog.txt")
            shutil.copy(source_path, dest_path)
            
            #copy config file from first subset #0, supposing the file is similar for the other subsets results
            #WARNING: todo this direct copy could lead to pb...
            source_path = os.path.join(subpathsList[0], "config.txt")
            dest_path = os.path.join(outputPath, "config.txt")
            shutil.copy(source_path, dest_path)
            
            #copy templates folder from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "templates")
            dest_path = os.path.join(outputPath, "templates")
            shutil.copytree(source_path, dest_path)
            
            #copy templates folder from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "templates_nocontinuum")
            dest_path = os.path.join(outputPath, "templates_nocontinuum")
            shutil.copytree(source_path, dest_path)
            
            #copy intermediate results if any
            self.copyIntermediateResultsDirs(subpathsList, outputPath, spcfileName="input.spectrumlist")
            
        except:
            print("WARNING: unable to merge some output files. Skipping that operation !")
            raw_input("\n\nWARNING: unable to merge some output files. will skip that operation.\nPress any key to continue...".format())
        
        
    def recombineAll(self):
        for f in self.flist:
            print("INFO: recombining from subsets: {}".format(f))
            self.recombine(f)
            
    def mergeCsvFiles(self, subpathsList, outPath, fileName):
        """
        input:  - subpathsList is a list of directories containing the 'fileName' file
                - outPath is a directory where the merged 'fileName' csv will be saved 
        """
        allLines = []
        
        for dirpath in subpathsList:
            fpath = os.path.join(dirpath, fileName)
            f = open(fpath, 'r')
            for line in f:
                lineStr = line.strip()
                if not lineStr.startswith('#'):
                    allLines.append(lineStr)
            f.close()
            
        fFilePath = os.path.join(outPath, fileName)
        fout = open(fFilePath, 'w')
        for l in allLines:
            fout.write(l)
            fout.write('\n')
        fout.close()
          
    def copyIntermediateResultsDirs(self, subpathsList, outPath, spcfileName):
        """
        """
        for dirpath in subpathsList:
            fpath = os.path.join(dirpath, spcfileName)
            f = open(fpath, 'r')
            for line in f:
                lineStr = line.strip()
                if not lineStr.startswith('#'):
                    spcName = lineStr.split("\t")[0]
                    spcNameNoExt = os.path.splitext(spcName)[0]
                    source_path = os.path.join(dirpath, spcNameNoExt)
                    dest_path = os.path.join(outPath, spcNameNoExt)
                    try:
                        shutil.copytree(source_path, dest_path)
                    except:
                        print("ERROR: unable to copy intermediate files for: {}".format(spcNameNoExt))
            f.close()
        
        
def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    

    parser.add_argument("-p", "--path", dest="path", default="./process-output/output_subsets",
                    help="path to the subset output directory to be combined")
    parser.add_argument("-o", "--outputpath", dest="outpath", default="./process-output/",
                    help="path to the destination output directory to fill with merged data") 
    
                    
    options = parser.parse_args()
    #print(options)

    if os.path.exists(options.path) :
        rpath = os.path.abspath(options.path)
        print('INFO: using full path: {}'.format(rpath))
        routpath = os.path.abspath(options.outpath)
        print('INFO: using full output path: {}'.format(routpath))
        precomb = processRecombine(path=rpath, outputPath=routpath)
    else :
        print("Error: invalid arguments")
        exit()
    
    
    
def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print(".")
    Main( sys.argv )
    

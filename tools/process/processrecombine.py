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
        
    def recombine(self, infoFilePath, enableSkipUnfinished=True):
        """
        recombine from the subset described by infoFilePath
        """
        subpathsList = []
        #load info file
        f = open(infoFilePath, 'r')
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                fullPath = os.path.join(self.outputPath, lineStr)
                subpathsList.append(fullPath)
            
        #create this config merged output path
        name = os.path.split(infoFilePath)[1]
        name_noext = os.path.splitext(name)[0]
        datasetOutputPath = os.path.join(self.outputPath, name_noext)
        if not os.path.exists(datasetOutputPath):
            os.mkdir(datasetOutputPath)
            
        #merge redshift.csv files
        self.mergeCsvFiles(subpathsList, datasetOutputPath, fileName="redshift.csv", verbose=True)
        print("WARNING: redshift.csv has been recombined!")
        
        #merge log files
        self.mergeCsvFiles(subpathsList, datasetOutputPath, fileName="log.txt")
        print("WARNING: log.txt has been recombined!")

        try:        
            #merge input spectrumlist files
            self.mergeCsvFiles(subpathsList, datasetOutputPath, fileName="input.spectrumlist")         
            
            #copy version file from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "versions.txt")
            dest_path = os.path.join(datasetOutputPath, "versions.txt")
            shutil.copy(source_path, dest_path)
            
            #copy json file from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "parameters.json")
            dest_path = os.path.join(datasetOutputPath, "parameters.json")
            shutil.copy(source_path, dest_path)
            
            #copy linecatalog file from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "linecatalog.txt")
            dest_path = os.path.join(datasetOutputPath, "linecatalog.txt")
            shutil.copy(source_path, dest_path)
            
            #copy config file from first subset #0, supposing the file is similar for the other subsets results
            #WARNING: todo this direct copy could lead to pb...
            source_path = os.path.join(subpathsList[0], "config.txt")
            dest_path = os.path.join(datasetOutputPath, "config.txt")
            shutil.copy(source_path, dest_path)
            
            #copy templates folder from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "templates")
            dest_path = os.path.join(datasetOutputPath, "templates")
            shutil.copytree(source_path, dest_path)
            
            #copy templates folder from first subset #0, supposing the file is similar for the other subsets results
            source_path = os.path.join(subpathsList[0], "templates_nocontinuum")
            dest_path = os.path.join(datasetOutputPath, "templates_nocontinuum")
            shutil.copytree(source_path, dest_path)
            
            #merge classification.csv files
            self.mergeCsvFiles(subpathsList, datasetOutputPath, fileName="classification.csv", verbose=True)
            print("WARNING: classification.csv has been recombined!")
            
            #merge stellar.csv files
            self.mergeCsvFiles(subpathsList, datasetOutputPath, fileName="stellar.csv", verbose=True)
            print("WARNING: stellar.csv has been recombined!")
  
        except Exception as e:
            if not enableSkipUnfinished:
                print("WARNING: unable to merge some output files. Skipping that operation !")
                input("\n\nWARNING: unable to merge some output files. will skip that operation.\nPress any key to continue...".format())
                raise e
            else:
                print(e)
                
                
        try:
            #copy intermediate results if any
            self.copyIntermediateResultsDirs(subpathsList, datasetOutputPath, spcfileName="input.spectrumlist")
                
        except Exception as e:
            if not enableSkipUnfinished:
                print("WARNING: unable to merge some INTERMEDIATE output files. Skipping that operation !")
                input("\n\nWARNING: unable to merge some output files. will skip that operation.\nPress any key to continue...".format())
                raise e
            else:
                print(e)
        
        #copy the cluster logs into the merged folder
        try:
            source_path = os.path.join(self.outputPath, "cluster_logs") #warning hardcoded, but not critical to avoid this...
            dest_path = os.path.join(datasetOutputPath, "cluster_logs")
            shutil.copytree(source_path, dest_path)
        except:
            print("INFO: skipped cluster_logs directory while recombining. Folder not found")

        #estimate per spectrum processing time stats
        try:
            clusterLogsDirPath= os.path.join(self.outputPath, "cluster_logs") #warning hardcoded, but not critical to avoid this...
            self.estimatePerSpectrumProcTime(clusterLogsDirPath, datasetOutputPath)
        except:
            print("INFO/Warning: skipped per spectrum processing time stats")

        
    def recombineAll(self):
        for f in self.flist:
            print("INFO: recombining from subsets: {}".format(f))
            self.recombine(f)
            
    def mergeCsvFiles(self, subpathsList, outPath, fileName, verbose=False):
        """
        input:  - subpathsList is a list of directories containing the 'fileName' file
                - outPath is a directory where the merged 'fileName' csv will be saved 
        """
        allLines = []
        headerCopied = False
        for dirpath in subpathsList:
            fpath = os.path.join(dirpath, fileName)
            f = open(fpath, 'r')
            if f==0:
                print("ERROR: unable to open file = {}".format(fpath))
            n=0
            for line in f:
                lineStr = line.strip()
                if not lineStr.startswith('#'):
                    allLines.append(lineStr)
                    n+=1
                elif not headerCopied:
                    headerCopied = True
                    allLines.append(lineStr)
            if verbose:
                print("INFO-merge: f={}, found lines n={}".format(fpath, n))            
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
                    procTag = lineStr.split("\t")[2]
                    source_path = os.path.join(dirpath, procTag)
                    dest_path = os.path.join(outPath, procTag)
                    try:
                        shutil.copytree(source_path, dest_path)
                    except:
                        print("ERROR: unable to copy intermediate files for: {}".format(procTag))
                        
            f.close()
            
    def estimatePerSpectrumProcTime(self, clusterLogDirPath, datasetOutputPath):
        """
        Look into the cluster log files and extract the processing time automatically
        current: uses some tags in the log files to find 'nspc' and 'proctime'. 
        current: (more): node identification is supposed to be the first line in the log file (from processhelper)
        """
        verbose = 1
        filOutPath = os.path.join(datasetOutputPath, "logProc.txt")
        fileOut = open(filOutPath, 'w')        
        fileOut.write("logSubset\tnodeName\tproctimesec\tnspc\tperSpcTime\n".format())
        
        #loop on the logOur files
        for i,file in enumerate(sorted(os.listdir(clusterLogDirPath))):
            if file.startswith("logOut"):  
                fullPath = os.path.join(clusterLogDirPath, clusterLogDirPath, file)
                if verbose:
                    print("INFO: estimating proc. time from: {}".format(file))
                # read nodeID, nSpc, processing time
                flog = open(fullPath)
                nodeStr = None 
                proctimeSeconds = -1  
                nSpc = 0
                for line in flog:
                    if nodeStr==None:                    
                        nodeStr=line.rstrip()
                    if "<proc-spc>" in line:
                        nSpc+=1
                    if "<proctime-seconds>" in line:
                        lineSplit = line.split("<proctime-seconds><")
                        if len(lineSplit)>1:
                            proctimeSeconds=int(lineSplit[1].split(">")[0])
                fileOut.write("{}\t{}\t{}\t{}\t{:.1f}\n".format(str(file), nodeStr, proctimeSeconds, nSpc, float(proctimeSeconds)/float(nSpc)))
        fileOut.close()
        
        return
        
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
    

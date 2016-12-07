# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 18:55:20 2016

@author: aschmitt
"""

## python 3 compatibility imports ###
from __future__ import division, print_function
from future_builtins import *
from builtins import range
from future.builtins.disabled import apply, cmp, file, raw_input, xrange
## ##### # ############# ####### ###

import os
import sys
import time
import argparse
import shutil
import json


class Report(object):
    def __init__(self, statsDirPath, outputDirPath, globalSuccessRateValue, parametersFilePath, configPath):
        self.statsdirpath = statsDirPath
        self.outputdirpath = outputDirPath
        self.parameterspath = parametersFilePath
        self.configpath = configPath
        self.templateDirPath = "/home/aschmitt/gitlab/cpf-redshift/tools/stats/latex_empty3"
        self.texFilePath = os.path.join(self.outputdirpath, "report_template.tex")
        
        #init report directory with template
        self.makeOutputDirFromTemplate(self.outputdirpath)

        #define the document details
        self.doc_date = time.strftime("%d/%m/%Y")  
        self.setTexPath("tpl_tobereplaced_formatteddate", self.doc_date) 

        
        #define figures relative paths
        self.perfMatrixSummaryPath = os.path.join(self.statsdirpath, 'performances', 'performance.png')
        stats_subset_path = 'stats_magmin-1.0magmax50.0_zmin-1.0zmax20.0_sfrmin-1.0sfrmax10000.0'
        self.successrateGlobalPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'stats_hist.png')
        self.scatterplotGlobalPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'filteredset_relzerr.png')
        self.scatterplotzcamczrefGlobalPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'filteredset_zcalc-vs-zref.png')
        self.scatterplotzoom0Path = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'filteredset_relzerr_zoom_0.png')
        self.scatterplotzoom1Path = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'filteredset_relzerr_zoom_1.png')
        self.scatterplotzoom2Path = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'filteredset_relzerr_zoom_2.png')
        self.histVersusMagPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'stats_versusMag_hist.png')
        self.histVersusRedshiftPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'stats_versusRedshift_hist.png')
        self.histVersusSfrPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'stats_versusSFR_hist.png')
        self.histVersusLOGFHalphaPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'stats_versusLogFHAlpha_hist.png')
        self.histVersusEbmvPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'stats_versusEBMV_hist.png')
        self.histVersusSigmaPath = os.path.join(self.statsdirpath, "{}".format(stats_subset_path), 'stats_versusSigma_hist.png')
        
        #define the global success rate figure caption
        self.globalSuccessRateBinCenter = "{:.2e}".format(globalSuccessRateValue[0])
        self.globalSuccessRateValue = "{:.2f}".format(globalSuccessRateValue[1])
        self.setTexPath("tpl_tobereplaced_globalsuccessrateBinCenter", self.globalSuccessRateBinCenter) 
        self.setTexPath("tpl_tobereplaced_globalsuccessrateValue", self.globalSuccessRateValue)  

        self.setTexPath("tpl_tobereplaced_fig_performancesummary", self.perfMatrixSummaryPath)  
        self.setTexPath("tpl_tobereplaced_fig_successrateglobal", self.successrateGlobalPath)
        self.setTexPath("tpl_tobereplaced_fig_scatterplotglobal", self.scatterplotGlobalPath)
        self.setTexPath("tpl_tobereplaced_fig_scatterplotzcalczrefglobal", self.scatterplotzcamczrefGlobalPath)
        self.setTexPath("tpl_tobereplaced_fig_scatterplotzoom0", self.scatterplotzoom0Path)
        self.setTexPath("tpl_tobereplaced_fig_scatterplotzoom1", self.scatterplotzoom1Path)
        self.setTexPath("tpl_tobereplaced_fig_scatterplotzoom2", self.scatterplotzoom2Path)
        self.setTexPath("tpl_tobereplaced_fig_histversus_mag", self.histVersusMagPath)
        self.setTexPath("tpl_tobereplaced_fig_histversus_redshift", self.histVersusRedshiftPath)
        self.setTexPath("tpl_tobereplaced_fig_histversus_sfr", self.histVersusSfrPath)
        self.setTexPath("tpl_tobereplaced_fig_histversus_logfhalpha", self.histVersusLOGFHalphaPath)
        self.setTexPath("tpl_tobereplaced_fig_histversus_ebmv", self.histVersusEbmvPath)
        self.setTexPath("tpl_tobereplaced_fig_histversus_sigma", self.histVersusSigmaPath)
        
        #fill the configuration table
        self.spectrumlist = os.path.split(self.getConfigVal('input'))[1]
        self.setTexPath("tpl_tobereplaced_param_spclist", self.spectrumlist, stype="text")
        self.spcPath = os.path.split(self.getConfigVal('spectrumdir'))[1]
        self.setTexPath("tpl_tobereplaced_param_spcdir", self.spcPath, stype="text")
        self.methodParam = self.getParameterVal('method')
        self.setTexPath("tpl_tobereplaced_param_method", self.methodParam, stype="text")
        self.spcCount = "-"
        self.setTexPath("tpl_tobereplaced_param_spccount", self.spcCount, stype="text")
        self.continuumRemovalMethod = self.getParameterVal('continuumRemoval', 'method')
        self.setTexPath("tpl_tobereplaced_param_continuummethod", self.continuumRemovalMethod, stype="text")
        
        #params path for the annexe...
        self.setTexPath("tpl_tobereplaced_jsonparamsfullpath", self.parameterspath)

        
        self.makePdf()
        
    def getParametersJsonString(self):
        filename = self.parameterspath
        f = open(filename, 'r')
        contents = f.read()
        return contents
    
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
    
    def makeOutputDirFromTemplate(self, outputDir):
        if os.path.exists( outputDir ) == False :        
            print("makedir: Output dir: "+outputDir)
            self.getTemplateFiles(outputDir)
        else:
            print("ERROR: Output dir: {} already exists. Please delete it manually and retry !\n\n".format(outputDir))
            exit(0)
        print("INFO: Using outputdir : {}".format(outputDir))
        
    def getTemplateFiles(self, destinationDirPath):
        _fromPath = self.templateDirPath
        _toPath = destinationDirPath
        shutil.copytree(_fromPath, _toPath)
        
    def setTexPath(self, tplString, destString, stype="path"):
        if not stype=="path":
            destString = destString.replace("_", "\_")
        
        _path = self.texFilePath
        f = open(_path,'r')
        contents = f.read()
        replaced_contents = contents.replace("{}".format(tplString), "{}".format(destString))
        f.close()
        f = open(_path,'w')
        f.write(replaced_contents)
        f.close()
        
    def makePdf(self):
        
        cwd = os.getcwd()
        os.chdir(self.outputdirpath)
        bashCommand = "pdflatex  --shell-escape report_template.tex".format()
        os.system(bashCommand) 
        os.chdir(cwd)
     

def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-s", "--statsDirPath", dest="statsDirPath", default="",
                    help="path to the stats directory")
    parser.add_argument("-o", "--outputDirPath", dest="outputDirPath", default="./report",
                    help="path to the output directory")


    options = parser.parse_args()

    #print(options)

    if os.path.exists(options.statsDirPath) :
        print('using full path: {0}'.format(options.statsDirPath))
        r = Report(statsDirPath=options.statsDirPath, outputDirPath=options.outputDirPath)
        
    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print("Report")
    Main( sys.argv )
           
   
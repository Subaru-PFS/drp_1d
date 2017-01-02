# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 08:25:11 2016

@author: aschmitt
"""

import os
import sys
import time
import argparse


#cmd_subfolder = os.path.abspath("/home/aschmitt/gitlab/cpf-redshift/tools/aview/") #use locally
#cmd_subfolder = os.path.abspath("/home/aschmitt/amazed_cluster/gitlab/cpf-redshift/tools/aview/") #use on hermes
#if cmd_subfolder not in sys.path:
#    print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
#    sys.path.insert(0, cmd_subfolder)
#import spectrum
import spectrumlist



class processHelper(object):
    def __init__(self, confpath, binpath, rootoutputpath, dividecount):
        self.logTagStr = "processHelper"
        self.configPath = confpath
        self.binPath = binpath
        self.baseoutputpath = rootoutputpath

        self.proc_date = time.strftime("%Y%m%d")          

        #self.opt_bracketing = "" #no bracketing
        self.opt_bracketing = "method" 
        self.bracketing_templatesRootPath = "/home/aschmitt/amazed_cluster/calibration/templates/"
        
        self.loadConfig()
        
        #prepare the working dir
        self.work_process_dir = os.path.abspath("work_process")
        if not os.path.exists(self.work_process_dir):
            os.mkdir(self.work_process_dir)        
        
        if dividecount>0 :
            outputPath = os.path.join(self.work_process_dir, "spectrumlist_subs_{}".format(dividecount))
            if not os.path.exists(outputPath):
                os.mkdir(outputPath)
            print('INFO: splitting using full path: {}'.format(outputPath))
            spclist = spectrumlist.Spectrumlist(self.config_spclistPath)
            self.subspclists = spclist.splitIntoSubsets(int(dividecount), outputPath)
            self.subrecombine_info = {}
            self.baseoutputpath = os.path.join(self.baseoutputpath, "output_subsets")
            if not os.path.exists(self.baseoutputpath):
                os.mkdir(self.baseoutputpath)
        else:
            self.subspclists=[]
            
                              
    def getConfigVal(self, tag):
        """
        returns the string value of the field corresponding to the tag in the config file 
        """
        strVal = "not found"
        filename = self.configPath
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
                        
    def loadConfig(self):
        """
        Load only the relevant stuff from the config file
        """
        print("INFO: Loading config file...")
        self.config_spclistPath = self.getConfigVal("input")
        print("INFO: input is : {}".format(self.config_spclistPath))
        self.config_spcdir = self.getConfigVal("spectrumdir")
        print("INFO: spectrum dir. is : {}".format(self.config_spcdir))
        self.config_method = self.getConfigVal("method")
        print("INFO: method is : {}".format(self.config_method))
        self.config_templatedir = self.getConfigVal("templatedir")
        print("INFO: templatedir is : {}".format(self.config_templatedir))
        self.config_linecatalog = self.getConfigVal("linecatalog")
        print("INFO: linecatalog is : {}".format(self.config_linecatalog))
        self.config_linecatalogconvert = "not found" not in self.getConfigVal("linecatalog-convert")
        print("INFO: linecatalogconvert is : {}".format(self.config_linecatalogconvert))
        
        self.config_parametersPath = self.getConfigVal("parameters")
        print("INFO: parameters path is : {}".format(self.config_parametersPath))

    def prepareBracketing(self):
        """
        prepare a dictionnary (keys=method, parameters, templates) of configurations to be tested
        """
        if self.opt_bracketing=="method":
            bracketing_method=["linemodel", "linemodel", "amazed0_3", "amazed0_1", "chisquare2solve", "chisquare2solve"]
            bracketing_parameterFilePath=["lm_rules.json", "lm_tplshape.json", "amazed0_3.json", "amazed0_1.json", "chi2_nc.json", "chi2_raw.json"]
            bracketing_templateCtlNames=["ExtendedTemplatesMarch2016_v2", "ExtendedTemplatesMarch2016_v2", "ContinuumTemplates_simulm2016", "ExtendedTemplatesMarch2016_v2", "ExtendedTemplatesMarch2016_v2", "ExtendedTemplatesMarch2016_v2"]   
            bracketing_templateCtlPath = []            
            for k in range(len(bracketing_templateCtlNames)):
                bracketing_templateCtlPath.append(os.path.join(self.bracketing_templatesRootPath, bracketing_templateCtlNames[k]))
                
        else:
            #toto implement cont. method bracketing or other params...
            pass
        
        dico = {"method":bracketing_method, "parameters":bracketing_parameterFilePath, "templates":bracketing_templateCtlPath}

        #print bracketing 
        nproc = len(dico["method"])     
        print("{:<12}{:<20}{:<30}{:<45}".format("", "Method", "Parameters", "Templates"))
        for k  in range(nproc):
            print("{:<12}{:<20}{:<30}{:<45}".format(k, dico["method"][k], dico["parameters"][k], dico["templates"][k]))
        print("\n")
        
        return dico
        
    def prepareArgumentString(self, overrideMethod="", overrideParamFile="", overrideTemplatesCtlg="", overrideSpclist="", overrideSpclistIndex=-1):
        """
        prepares the args for the amazed bin call
        - default args come from the input config file loaded previously
        - overrided args allow for bracketing and use of divided spclists
        """
        argStr = ""
        if overrideSpclist=="":
            splListPath = self.config_spclistPath
        else:
            splListPath = overrideSpclist
        argStr = "{} --input {}".format(argStr, splListPath)
        argStr = "{} --spectrumdir {}".format(argStr, self.config_spcdir)

        if overrideMethod=="":
            argMethod = self.config_method
        else:
            argMethod = overrideMethod
        argStr = "{} --method {}".format(argStr, argMethod)
            
        if overrideTemplatesCtlg=="":
            argTpl = self.config_templatedir
        else:
            argTpl = overrideTemplatesCtlg
        argStr = "{} --templatedir {}".format(argStr, argTpl)
        
        argStr = "{} --linecatalog {}".format(argStr, self.config_linecatalog)
        if self.config_linecatalogconvert:
            argStr = "{} --linecatalog-convert".format(argStr)

        if overrideParamFile=="":
            argStr = "{} --parameters {}".format(argStr, self.config_parametersPath)
            argParamsLabelWExt = os.path.split(self.config_parametersPath)[1]
            argParamsLabel = os.path.splitext(argParamsLabelWExt)[0]
        else:
            #Note: looks for the bracketing param files in the original folder from config file
            overrideParamFolder = os.path.split(self.config_parametersPath)[0]
            overrideParamFileFullPath = os.path.join(overrideParamFolder, overrideParamFile)
            argStr = "{} --parameters {}".format(argStr, overrideParamFileFullPath)
            argParamsLabel = os.path.splitext(overrideParamFile)[0]            
            
        #warning, hardcoded: always use only 1 proc. thread        
        argStr = "{} --thread-count {}".format(argStr, "1")
        
           
        outputName_base = "amazed_{}_m-{}_p-{}".format(self.proc_date, argMethod, argParamsLabel)
        if overrideSpclistIndex==-1:
            outputName = outputName_base
        else:
            argSubIndex = "_sub-{}".format(overrideSpclistIndex)
            outputName = "{}{}".format(outputName_base, argSubIndex)
    
        outputPath = os.path.join(self.baseoutputpath, outputName)
        argStr = "{} --output {}".format(argStr, outputPath)
        
        if not overrideSpclistIndex==-1:
            if not outputName_base in self.subrecombine_info.keys():
                self.subrecombine_info[outputName_base] = []
            self.subrecombine_info[outputName_base].append(outputPath)
            
        return argStr        
        
    def exportSubRecombineInfo(self):
        
        for k, key in enumerate(self.subrecombine_info.keys()):
            print("INFO: exporting recombine info for {}".format(key))
            outputName = "{}.info".format(key)
            outputPath = os.path.join(self.baseoutputpath, outputName)
            fout = open(outputPath, 'w')
            for ksub, outPathSub in enumerate(self.subrecombine_info[key]):            
                #print("INFO: \toutPathSub = {}".format(outPathSub))
                fout.write("{}".format(outPathSub))
                fout.write("\n")
            fout.close()
            

    def doProcessLocal(self):
        """
        process on the local machine
        """
        print("doProcessLocal...")
        print("WARNING: LOCAL processing is DEPRECATED! not up to date !!!!, aborting...")
        return        
        
        if self.opt_bracketing == "method":
            bracketing = self.prepareBracketing()        
        else:
            bracketing = ""
                        
        for k in range(len(bracketing["method"])):  
            if bracketing=="":              
                argStr = self.prepareArgumentString()
            else:
                argStr = self.prepareArgumentString(overrideMethod = bracketing["method"][k], overrideParamFile = bracketing["parameters"][k])
            #print("INFO: arg string is {}".format(argStr))
            
            bashCommand = "{} {}".format(self.binPath, argStr)
            print("INFO: bashCommand is {}".format(bashCommand))
            os.system(bashCommand) 
    
    def doProcessQsub(self, dryrun=False):
        """
        prepare the .sh files to be sent to be qsub'd to the cluster
        1. prepare bracketing parameters
        2. loop on the bracketing index (method or other bracketing params)
        3. loop on the sub spclists
        4. submit pbs files to batch-cluster
        5. export the recombine info for post-processing recombination of the spclist-subsets
        """
        print("doProcessQsub...")
        
        
        if self.opt_bracketing == "method":
            bracketing = self.prepareBracketing()    
            nproc = len(bracketing["method"])
        else:
            bracketing = ""
            nproc=1            
                                    
        for k in range(nproc):  
            
            if len(self.subspclists)==0:
                nsubspclist=1
            else:
                nsubspclist=len(self.subspclists)
                
            for ksubs in range(nsubspclist):
            
                if len(self.subspclists)==0:
                    overrideSpclist=""
                else:
                    overrideSpclist=self.subspclists[ksubs]
                    
                if bracketing=="":              
                    argStr = self.prepareArgumentString(overrideSpclist=overrideSpclist)
                else:
                    argStr = self.prepareArgumentString(overrideMethod = bracketing["method"][k], 
                                                        overrideParamFile = bracketing["parameters"][k], 
                                                        overrideTemplatesCtlg = bracketing["templates"][k],
                                                        overrideSpclist = overrideSpclist,
                                                        overrideSpclistIndex = ksubs)
                #print("INFO: arg string is {}".format(argStr))
                
                fpbsname = "pbs_{}_{}.sh".format(k, ksubs)
                fpbspath = os.path.join(self.work_process_dir, fpbsname)
                f = open(fpbspath, 'w')
                #f.write("#PBS -t	1")
                #f.write("\n")
                f.write("#PBS -N amazed_{}_{}".format(k, ksubs))
                f.write("\n")
                f.write("#PBS -l nodes=1:ppn=1")
                f.write("\n")
                f.write("#PBS -l walltime=319:00:00")
                f.write("\n")
                f.write("\n")
                f.write("cd $PBS_O_WORKDIR")
                f.write("\n")
                #f.write("taskID=$PBS_ARRAYID")
                #f.write("\n")
                f.write("\n")
                f.write("logErr={}".format(os.path.join(self.baseoutputpath, "logErr_{}-{}.txt".format(k, ksubs) )))
                f.write("\n")
                f.write("logOut={}".format(os.path.join(self.baseoutputpath, "logOut_{}-{}.txt".format(k, ksubs) )))
                f.write("\n")
                f.write("\n")
                
                argStr = "{} 2>${{logErr}}>${{logOut}}".format(argStr)
    
                f.write("{} {}".format(self.binPath, argStr))
                f.write("\n")
                
                #f.write("")
                f.close()            
                       
                if not dryrun:
                    bashCommand = "qsub {}".format(fpbspath)
                    print("INFO: bashCommand is {}\n\n".format(bashCommand))
                    os.system(bashCommand) 
                    
        if not len(self.subspclists)==0: 
            self.exportSubRecombineInfo()
        
        
def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    

    parser.add_argument("-b", "--binPath", dest="binPath", default="/home/aschmitt/amazed_cluster/gitlab/amazed/bin/amazed-0.0.0",
                    help="path to the amazed binary file") 
    parser.add_argument("-c", "--configPath", dest="configPath", default="",
                    help="path to the config file to be loaded") 
                    
    #option to split the spectrumlist
    parser.add_argument("-d", "--divide", dest="divide", action='store_true',
                    help="enable spectrumlist dividing/splitting into subsets")
    parser.add_argument("-n", "--dividecount", dest="dividecount", default=100,
                    help="dividing/splitting count") 
    parser.add_argument("-o", "--outputPath", dest="outputPath", default="output_process",
                    help="path to the root output path") 
                    
    #option to do it local or qsub
    parser.add_argument("-l", "--local", dest="local", action='store_true',
                    help="enable process locally, default (False) is to send the jobs to the cluster through qsub")
                    
                    
                    
    options = parser.parse_args()
    #print(options)

    if os.path.exists(options.configPath) :
        rpath = os.path.abspath(options.configPath)
        print('INFO: using config full path: {}'.format(rpath))
        rbinpath = os.path.abspath(options.binPath)

        routputpath = os.path.abspath(options.outputPath)
        if not os.path.exists(routputpath):
            os.mkdir(routputpath)
        
        if options.divide:
            dividecount = options.dividecount
        else:
            dividecount =-1
            
        cp = processHelper(confpath=rpath, binpath=rbinpath, rootoutputpath=routputpath, dividecount=dividecount)
            
        if not options.local:
            cp.doProcessQsub()
        else:
            cp.doProcessLocal()
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
    

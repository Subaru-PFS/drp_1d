# -*- coding: utf-8 -*-
"""
Created on Sat Dec 24 08:25:11 2016

@author: aschmitt
"""

import os
import sys
import time
from datetime import datetime
import argparse

import math
import processParamCheck
import json

try:
    #cmd_subfolder = os.path.abspath("/home/aschmitt/gitlab/cpf-redshift/tools/aview/") #use locally
    cmd_subfolder = os.path.abspath("/home/aschmitt/amazed_cluster/gitlab/cpf-redshift/tools/aview_hardcopy/") #use on hermes
    if cmd_subfolder not in sys.path:
        print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
        sys.path.insert(0, cmd_subfolder)
        
    import reference
except ImportError:
    print("Import ERROR: unable to load the reference package. some functionnalities will be unavailable...")

import spectrumlist




class processHelper(object):
    def __init__(self, confpath, binpath, rootoutputpath, dividecount, opt_bracketing, bracketing_templatesRootPath, refpath):
        self.logTagStr = "processHelper"
        self.ready = False
        self.configPath = confpath
        self.binPath = binpath
        self.baseoutputpath = rootoutputpath
        
        self.logsPath = os.path.join(self.baseoutputpath, "cluster_logs")
        if not os.path.exists(self.logsPath):
            os.mkdir(self.logsPath)
            
        self.refPath = refpath
        if self.refPath == "":
            self.enableProcessAtZ = 0
        else:
            self.enableProcessAtZ = 1
        if self.enableProcessAtZ and not dividecount==1:
            print("ERROR: incompatible options : enableProcessAtZ and dividecount!=1. Aborting.")
            return
        if self.enableProcessAtZ:
            self.refcatalog = reference.Reference(referencepath=self.refPath, rtype="simple")

        self.proc_date = time.strftime("%Y%m%d")          

        self.opt_bracketing = opt_bracketing #choices = "", "method"
        self.bracketing_templatesRootPath = bracketing_templatesRootPath
        
        ret = self.loadConfig()
        if not ret:
            print("ERROR: load config failed.")
            return
        
        #prepare the working dir
        fn_datestr = datetime.now().strftime('%Y%m%d_%H%M%S')
        self.work_process_dir = os.path.abspath("process-work_{}".format(fn_datestr))
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
            self.subsetsRelPath = "output_subsets"
            self.baseoutputpath = os.path.join(self.baseoutputpath, self.subsetsRelPath)
            if not os.path.exists(self.baseoutputpath):
                os.mkdir(self.baseoutputpath)
        else:
            self.subspclists=[]
            
        self.ready = True
            
                              
    def getConfigVal(self, tag):
        """
        returns the string value of the field corresponding to the tag in the config file 
        """
        verbose=0
        strVal = "not found"
        filename = self.configPath
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            if lineStr.startswith("#"):
                continue
            if verbose:
                print("DEBUG: lineStr = {}".format(lineStr))
            data = lineStr.split("=")
            if verbose:
                print("DEBUG: data splitted = {}".format(data))            
                print("DEBUG: len(data) = {}".format(len(data)))
            if(len(data) >=2):
                if verbose:
                    print("DEBUG: tag={}".format(tag))
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
        if not os.path.exists(self.config_spclistPath):
            print("ERROR: spclist file does not exist! Aborting...")
            return False
            
        self.config_spcdir = self.getConfigVal("spectrumdir")
        print("INFO: spectrum dir. is : {}".format(self.config_spcdir))
        if not os.path.exists(self.config_spcdir):
            print("ERROR: spcdir file does not exist! Aborting...")
            return False
            
        self.config_method = self.getConfigVal("method")
        print("INFO: method is : {}".format(self.config_method))
        self.config_templatedir = self.getConfigVal("templatedir")
        print("INFO: templatedir is : {}".format(self.config_templatedir))
        if not os.path.exists(self.config_templatedir):
            print("ERROR: templatedir file does not exist! Aborting...")
            return False
        
        self.config_calibrationdir = self.getConfigVal("calibrationdir")
        print("INFO: calibrationdir is : {}".format(self.config_calibrationdir))
        if not os.path.exists(self.config_calibrationdir):
            print("ERROR: calibration dir does not exist! Aborting...")
            return False
            
                
        self.config_zclassifierdir = self.getConfigVal("zclassifierdir")
        print("INFO: reliabilitydir is : {}".format(self.config_zclassifierdir))
        enableSkipReliability = True
        if not os.path.exists(self.config_zclassifierdir):
            if not enableSkipReliability:
                print("ERROR: zclassifierdir dir does not exist! Aborting...")
                return False
            else:
                print("ERROR: zclassifierdir dir does not exist! Skipping...")
            
        self.config_linecatalog = self.getConfigVal("linecatalog")
        print("INFO: linecatalog is : {}".format(self.config_linecatalog))
        if not os.path.exists(self.config_linecatalog):
            print("ERROR: linecatalog file does not exist! Aborting...")
            return False
            
        self.config_linecatalogconvert = "not found" not in self.getConfigVal("linecatalog-convert")
        print("INFO: linecatalogconvert is : {}".format(self.config_linecatalogconvert))
        
        self.config_parametersPath = self.getConfigVal("parameters")
        print("INFO: parameters path is : {}".format(self.config_parametersPath))
        if not os.path.exists(self.config_parametersPath):
            print("ERROR: parameters file does not exist! Aborting...")
            return False
        if 1:
            pcheck = processParamCheck.processParamCheck(self.config_parametersPath)
            retDeprecatedKeyword = pcheck.checkAll()
            if retDeprecatedKeyword:
                print("INFO: aborting...")
                return False
            
        self.config_saveintermediateresults = self.getConfigVal("saveintermediateresults")
        print("INFO: saveintermediateresults is : {}".format(self.config_saveintermediateresults))
        if not (self.config_saveintermediateresults=="no" or self.config_saveintermediateresults=="global" or self.config_saveintermediateresults=="default" or self.config_saveintermediateresults=="linemeas" or self.config_saveintermediateresults=="all"):
            print("ERROR: config_saveintermediateresults bad value (={})...".format(self.config_saveintermediateresults))
            return False
       
 
        self.config_verbose = self.getConfigVal("verbose")
        print("INFO: verbose is : {}".format(self.config_verbose)) 
        if not (self.config_verbose=="0" or self.config_verbose=="1" or self.config_verbose=="2"):
            print("ERROR: config_verbose bad value (={})...".format(self.config_verbose))
            return False
                
        self.config_linemeascatalogpath = self.getConfigVal("linemeascatalog")
        print("INFO: linemeascatalogpath is : {}".format(self.config_linemeascatalogpath))

       
        return True
           
        

    def prepareBracketing(self):
        """
        prepare a dictionnary (keys=method, parameters, templates) of configurations to be tested
        """
        if self.opt_bracketing=="method":
            bracketing_method=["linemodel", "linemodel", "amazed0_3", "amazed0_1", "chisquare2solve", "chisquare2solve"]
            bracketing_parameterFileName=["lm_rules.json", "lm_tplshape.json", "amazed0_3.json", "amazed0_1.json", "chi2_nc.json", "chi2_raw.json"]
            overrideParamFolder = os.path.split(self.config_parametersPath)[0]
            bracketing_parameterFilePath = []
            for k in range(len(bracketing_parameterFileName)):
                bracketing_parameterFilePath.append(os.path.join(overrideParamFolder, bracketing_parameterFileName[k]))
            bracketing_templateCtlNames=["ExtendedTemplatesMarch2016_v2", "ExtendedTemplatesMarch2016_v2", "ContinuumTemplates_simulm2016", "ExtendedTemplatesMarch2016_v2", "ExtendedTemplatesMarch2016_v2", "ExtendedTemplatesMarch2016_v2"]   
            bracketing_templateCtlPath = []            
            for k in range(len(bracketing_templateCtlNames)):
                bracketing_templateCtlPath.append(os.path.join(self.bracketing_templatesRootPath, bracketing_templateCtlNames[k]))
        elif self.opt_bracketing=="method-euclid":
            bracketing_method=["linemodel", 
                               "linemodel", 
                               "linemodel", 
                               "blindsolve", 
                               "chisquare2solve", 
                               "chisquare2solve"]
            bracketing_parameterFileName=["fullmodel_elv370alv530.json", 
                                          "fullmodel.json", 
                                          "linemodel_elv370alv530.json", 
                                          "blindsolve.json", 
                                          "chi2_nc.json", 
                                          "chi2_raw.json"]
            overrideParamFolder = os.path.split(self.config_parametersPath)[0]
            bracketing_parameterFilePath = []
            for k in range(len(bracketing_parameterFileName)):
                bracketing_parameterFilePath.append(os.path.join(overrideParamFolder, bracketing_parameterFileName[k]))
            bracketing_templateCtlNames=["ContinuumTemplates_simulm2016Extended_dustfree201702", 
                                         "ContinuumTemplates_simulm2016Extended_dustfree201702", 
                                         "ContinuumTemplates_simulm2016Extended_dustfree201702", 
                                         "ExtendedTemplatesJan2017_v3", 
                                         "ExtendedTemplatesJan2017_v3", 
                                         "ExtendedTemplatesJan2017_v3"]   
            bracketing_templateCtlPath = []            
            for k in range(len(bracketing_templateCtlNames)):
                bracketing_templateCtlPath.append(os.path.join(self.bracketing_templatesRootPath, bracketing_templateCtlNames[k]))
        elif self.opt_bracketing=="tplcat_analysis-pfs":
            bracketing_method=["chisquare2solve", "chisquare2solve", "chisquare2solve", "chisquare2solve", "chisquare2solve", "chisquare2solve"]
            bracketing_parameterFileName=["chisquare_rest_lbda4z0.json", "chisquare_rest_lbda4z0p5.json", "chisquare_rest_lbda4z1.json", "chisquare_rest_lbda4z1p5.json", "chisquare_rest_lbda4z2.json", "chisquare_rest_lbda4z2p5.json"]
            overrideParamFolder = os.path.split(self.config_parametersPath)[0]
            bracketing_parameterFilePath = []
            for k in range(len(bracketing_parameterFileName)):
                bracketing_parameterFilePath.append(os.path.join(overrideParamFolder, bracketing_parameterFileName[k]))
                
            #bracketing_templateCtlNames=["BC03_sdss_tremonti21", "BC03_sdss_tremonti21", "BC03_sdss_tremonti21", "BC03_sdss_tremonti21", "BC03_sdss_tremonti21", "BC03_sdss_tremonti21"]
            bracketing_templateCtlNames=["legac_z0_697", "legac_z0_697", "legac_z0_697", "legac_z0_697", "legac_z0_697", "legac_z0_697"]   
            bracketing_templateCtlPath = []            
            for k in range(len(bracketing_templateCtlNames)):
                bracketing_templateCtlPath.append(os.path.join(self.bracketing_templatesRootPath, bracketing_templateCtlNames[k]))
        elif self.opt_bracketing=="priors_optimization-euclid":
            
            #rawParamFolder = os.path.split(self.config_parametersPath)[0]
            tmp_subfolder = "tmp_parameters_brackets"
            tmp_bracketing_parameters_dir = os.path.join(self.work_process_dir, tmp_subfolder)
            if not os.path.exists(tmp_bracketing_parameters_dir):
                os.mkdir(tmp_bracketing_parameters_dir)
            
            nhaempriorSs = [-1, 1, 2, 3, 4, 5, 6]
            selpps = [-1, 1e-1, 1e-2, 1e-3, 1e-6, 1e-12]
            
            bracketing_parameterFilePath=[]
            for nhaempriorS in nhaempriorSs:
                for selpp in selpps:
                    argParamsLabelWExt = os.path.split(self.config_parametersPath)[1]
                    argParamsLabel = os.path.splitext(argParamsLabelWExt)[0]
                    
                    if nhaempriorS<0:
                        nhaempsStr="NO"
                    else:
                        nhaempsStr = ("{}".format(nhaempriorS)).replace("-", "m").replace(".", "p")
                    if selpp<0:
                        selppStr="NO"
                    else:
                        selppStr = ("{:.3e}".format(selpp)).replace("-", "m").replace(".", "p").replace("+", "")
                    overrideParamFile = "{}_nhaemps{}_selpp{}.json".format(argParamsLabel, nhaempsStr, selppStr)                  
                    overrideParamFileFullPath = os.path.join(tmp_bracketing_parameters_dir, overrideParamFile)
                    
                    fin=open(self.config_parametersPath)
                    data_json = json.load(fin)
                    #data_json['linemodelsolve']['linemodel']['stronglinesprior']=selpp
                    data_json['linemodelsolve']['linemodel']['haprior']=selpp
                    data_json['linemodelsolve']['linemodel']['euclidnhaemittersStrength']=nhaempriorS
                    fin.close()
                    # Writing JSON data
                    with open(overrideParamFileFullPath, 'w') as fout:
                        json.dump(data_json, fout)
                    fout.close()
                    
                                        
                    bracketing_parameterFilePath.append(overrideParamFileFullPath)
            
            n_brackets = len(bracketing_parameterFilePath)
            bracketing_method = [self.config_method for a in range(n_brackets)]
            bracketing_templateCtlPath = [self.config_templatedir for a in range(n_brackets)]   

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
        
    def prepareArgumentString(self, overrideMethod="", overrideParamFilePath="", overrideTemplatesCtlg="", overrideSpclist="", overrideSpclistIndex=-1):
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
        
        argCalib = self.config_calibrationdir
        argStr = "{} --calibrationdir {}".format(argStr, argCalib)

        argreliability = self.config_zclassifierdir        
        argStr = "{} --zclassifierdir {}".format(argStr, argreliability)
        
        argStr = "{} --linecatalog {}".format(argStr, self.config_linecatalog)
        if self.config_linecatalogconvert:
            argStr = "{} --linecatalog-convert".format(argStr)

        if overrideParamFilePath=="":
            argStr = "{} --parameters {}".format(argStr, self.config_parametersPath)
            argParamsLabelWExt = os.path.split(self.config_parametersPath)[1]
            argParamsLabel = os.path.splitext(argParamsLabelWExt)[0]
        else:         
            argStr = "{} --parameters {}".format(argStr, overrideParamFilePath)
            argParamsLabelWExt = os.path.split(overrideParamFilePath)[1]
            argParamsLabel = (os.path.splitext(argParamsLabelWExt)[0])  
            
        argStr = "{} --saveintermediateresults {}".format(argStr, self.config_saveintermediateresults) 
            
        #warning, hardcoded: always use only 1 proc. thread        
        argStr = "{} --thread-count {}".format(argStr, "1")
                
        argStr = "{} --verbose {}".format(argStr, self.config_verbose)
       
 
        if not os.path.exists(self.config_linemeascatalogpath) and (self.config_linemeascatalogpath!="not found"):
            print("ERROR: linemeascatalog bad value (={})...".format(self.config_linemeascatalogpath))
            return False
        elif (self.config_linemeascatalogpath!="not found"):
            argStr = "{} --linemeascatalog {}".format(argStr, self.config_linemeascatalogpath)




        if self.enableProcessAtZ:
            #Warning: hardcoded read first spc entry in the spclist, assuming it's the spc nametag 
            fspc = open(splListPath)
            for line in fspc:
                lineStr = line.strip()
                if not lineStr.startswith('#'):
                    data = lineStr.split("\t")
                    if len(data)<2:
                        data = lineStr.split(" ")
                    data = [r for r in data if r != '']
                    spcName = data[0]
                    break
            zprocess = self.refcatalog.getRedshift(spcName)
            zstep = 1e-6
            zrangehalf = 0.0
            #zrangehalf = 1e-3 #enable exploring the redshift range around zref
            argStr = """{} --redshiftrange "{} {} {}" """.format(argStr, zprocess-zrangehalf, zprocess+zrangehalf, zstep) #ex, format: redshiftrange=0.0 5.0 0.0001
           
        outputName_base = "amazed_{}_m-{}_p-{}".format(self.proc_date, argMethod, argParamsLabel)
        if overrideSpclistIndex==-1:
            outputName = outputName_base
        else:
            argSubIndex = "_sub-{}".format(overrideSpclistIndex)
            outputName = "{}{}".format(outputName_base, argSubIndex)
    
        outputPath = os.path.join(self.baseoutputpath, outputName)
        argStr = "{} --output {}".format(argStr, outputPath)
        
        if not overrideSpclistIndex==-1:
            relativePath = os.path.join(self.subsetsRelPath, outputName)
            if not outputName_base in self.subrecombine_info.keys():
                self.subrecombine_info[outputName_base] = []
            self.subrecombine_info[outputName_base].append(relativePath)
            
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
    
    def doProcess(self, dryrun=False, local="local"):
        """
        prepare the .sh files to be sent to be qsub'd to the cluster
        1. prepare bracketing parameters
        2. loop on the bracketing index (method or other bracketing params)
        3. loop on the sub spclists
        4. submit pbs files to batch-cluster
        5. export the recombine info for post-processing recombination of the spclist-subsets
        """
        print("doProcessQsub...")
        if not self.ready:
            print("processhelper not ready, aborting...")
            return
        
        if self.opt_bracketing == "method" or self.opt_bracketing == "method-euclid" or self.opt_bracketing == "tplcat_analysis-pfs" or self.opt_bracketing == "priors_optimization-euclid":
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
                
            skip_n_first_jobs = -1
            for ksubs in range(nsubspclist):
                if ksubs<skip_n_first_jobs:
                    continue
                
                if len(self.subspclists)==0:
                    overrideSpclist=""
                    overrideSpclistIndex = -1
                else:
                    overrideSpclist=self.subspclists[ksubs]
                    overrideSpclistIndex = ksubs
                    
                if bracketing=="":              
                    argStr = self.prepareArgumentString(overrideSpclist=overrideSpclist,
                                                        overrideSpclistIndex = overrideSpclistIndex)
                else:
                    argStr = self.prepareArgumentString(overrideMethod = bracketing["method"][k], 
                                                        overrideParamFilePath = bracketing["parameters"][k], 
                                                        overrideTemplatesCtlg = bracketing["templates"][k],
                                                        overrideSpclist = overrideSpclist,
                                                        overrideSpclistIndex = overrideSpclistIndex)
                #print("INFO: arg string is {}".format(argStr))
                
                #creating .sh file
                logepath = os.path.join(self.logsPath, "logErr_{}-{}.txt".format(k, ksubs) )
                logopath = os.path.join(self.logsPath, "logOut_{}-{}.txt".format(k, ksubs) )
                if local=="cluster_lam": #cluster PBS files
                    fpbsname = "pbs_{}_{}.sh".format(k, ksubs)
                    fpbspath = os.path.join(self.work_process_dir, fpbsname)
                    f = open(fpbspath, 'w')
                    #f.write("#PBS -t	1")
                    #f.write("\n")
                    f.write("#PBS -N amz_{}_{}".format(k, ksubs))
                    f.write("\n")
                    f.write("#PBS -l nodes=1:ppn=1")
                    f.write("\n")
                    f.write("#PBS -l mem=2GB") # request 2GB of memory
                    f.write("\n")
                    f.write("#PBS -l walltime=30:00:00")
                    f.write("\n")
                    f.write("\n")
                    f.write("cd $PBS_O_WORKDIR")
                    f.write("\n")
                    #f.write("taskID=$PBS_ARRAYID")
                    #f.write("\n")
                    f.write("\n")
                    f.write("logErr={}".format(logepath))
                    f.write("\n")
                    f.write("logOut={}".format(logopath))
                    f.write("\n")
                    f.write("\n")
                    f.write("hostname>${logOut}")
                    f.write("\n")
                    f.write("date>>${logOut}")
                    f.write("\n")
                    
                    f.write("\n")
                    argStr = "{} 2>${{logErr}}>>${{logOut}}".format(argStr)
        
                    f.write("{} {}".format(self.binPath, argStr))
                    f.write("\n")
                    
                    #f.write("")
                    f.close()  
                elif local=="cluster_in2p3":
                    fpbsname = "pbs_{}_{}.sh".format(k, ksubs)
                    fpbspath = os.path.join(self.work_process_dir, fpbsname)
                    f = open(fpbspath, 'w')
                    
                    f.write("\n")
        
                    f.write("cd /sps/euclid/Users/schmitt/amazed_cluster")
                    f.write("\n")
                    
                    f.write("{} {}".format(self.binPath, argStr))
                    f.write("\n")
                    
                    #f.write("")
                    f.close() 
                else: #local bash files
                    fshname = "amazed-proc_{}_{}.sh".format(k, ksubs)
                    fshpath = os.path.join(self.work_process_dir, fshname)
                    f = open(fshpath, 'w')
                    f.write("{} {}".format(self.binPath, argStr))
                    f.write("\n")
                    f.close() 
                    
                if not dryrun:
                    if local=="cluster_lam":
                        bashCommand = "qsub {}".format(fpbspath)
                    elif local=="cluster_in2p3":
                        options = "-l sps=1"
                        options += " " ##
                        options += "-l os=cl7"
                        options += " " ##
                        
                        
                        options += "-P P_euclid_spe"
                        options += " " ##
                        options += "-N amz_{}_{}".format(k, ksubs)
                        options += " " ##
                        #options += "-cwd {}".format()
                        options += " " ##
                        options += "-e {}".format(logepath)
                        options += " " ##
                        options += "-o {}".format(logopath)
                        options += " " ##
                        
                        bashCommand = "qsub {} {}".format(options, fpbspath)
                        
                    else:
                        bashCommand = "bash {}".format(fshpath)
                    print("INFO: bashCommand is {}\n\n".format(bashCommand))
                    os.system(bashCommand) 
                    
                    #try to sleep n seconds in order to avoid problems at submission
                    time.sleep(0.25)
                    
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
    parser.add_argument("-o", "--outputPath", dest="outputPath", default="process-output",
                    help="path to the root output path") 
                    
    #option to do it local or qsub
    parser.add_argument("-l", "--local", dest="local", default='cluster_lam',
                    help="dryrun (skip run), local (run locally), 'cluster_lam' or 'cluster_in2p3' to run on a cluster through qsub")
                 
    parser.add_argument("-m", "--methodBracketing", dest="methodBracketing", default="",
                help="bracketing of the parameters : choose between '' or 'method' or 'method-euclid' or 'tplcat_analysis-pfs' or 'priors_optimization-euclid' ") 
                
    parser.add_argument("-t", "--bracketingTemplatesRootPath", dest="brackTemplatesRootPath", default="/home/aschmitt/amazed_cluster/calibration/templates/",
                help="path to the templates root path used for the bracketing")
                
    parser.add_argument("-r", "--refPath", dest="refPath", default="",
                help="path to the reference file, used for processing at z for each source. Forces spclist splitting with count=1") 
    
                    
    options = parser.parse_args()
    #print(options)

    if os.path.exists(options.configPath) :
        rpath = os.path.abspath(options.configPath)
        print('INFO: using config full path: {}'.format(rpath))
        
        rbinpath = os.path.abspath(options.binPath)
        if not os.path.exists(rbinpath):
            print("ERROR: amazed binary file not found. Aborting")
            exit()
        else:
            print('INFO: using bin full path: {}'.format(rbinpath))
        
        routputpath = os.path.abspath(options.outputPath)
        if not os.path.exists(routputpath):
            os.mkdir(routputpath)
        print('INFO: using output full path: {}'.format(routputpath))
        
        rrefPath = ""
        if options.divide and options.refPath=="":
            dividecount = int(options.dividecount)
        elif not options.refPath=="":
            print("INFO: ref file path has been provided: spclist splitting forced with dividecount=1")
            dividecount = 1
            rrefPath = options.refPath
            if not os.path.exists(rrefPath):
                print("ERROR: reference file not found. Aborting")
                exit()
        else:
            dividecount = -1
                       
        opt_bracketing = options.methodBracketing
        print("INFO: cmdline arg bracketing = {}".format(opt_bracketing))
        
        bracketing_templatesRootPath = ""
        if not opt_bracketing=="" and not opt_bracketing=="priors_optimization-euclid":
            bracketing_templatesRootPath = os.path.abspath(options.brackTemplatesRootPath)
            if not os.path.exists(bracketing_templatesRootPath):
                print("ERROR: bracketing templates root path not found: {}\nAborting".format(bracketing_templatesRootPath))
                exit()
                
        cp = processHelper(confpath=rpath, binpath=rbinpath, rootoutputpath=routputpath, dividecount=dividecount, opt_bracketing=opt_bracketing, bracketing_templatesRootPath=bracketing_templatesRootPath, refpath=rrefPath)
            
        if options.local == "dryrun":
            dryrun = True
        else:
            dryrun=False
        cp.doProcess(dryrun=dryrun, local=options.local)
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
    

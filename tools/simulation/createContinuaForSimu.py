# -*- coding: utf-8 -*-
"""
"""
import time
import sys
import os
import shutil
import math
import random
import inspect
import numpy as np

subfolder = "../aview"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],subfolder)))
if cmd_subfolder not in sys.path:
    #print("inserting sys path : cmd_subfolder = {}".format(cmd_subfolder))
    sys.path.insert(0, cmd_subfolder)
import spectrum
import utilbins
import catalog as ctlg

def loadRef(fname):
    """
    load the ref file data
    """ 
    dataArray = []
    f = open(fname)
    for line in f:
        lineStr = line.strip()
        if not lineStr.startswith('#'):
            #print lineStr
            data = lineStr.split("\t")
            if len(data)<2:
                data = lineStr.split(" ")
            data = [r for r in data if r != '']
            #print len(data)
            if 1:
                if(len(data) >= 8): #SIMUEuclid2016
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = str(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = float(data[6])
                    d7 = float(data[7])
                    d = [d0, d1, d2, d3, d4, d5, d6, d7]
                    dataArray.append(d) 

            
    f.close()
    return dataArray
    
    

def Main( argv ) :	
    try:
        output_directory = "./simulation_output_tpl/".format()
        if os.path.exists(output_directory):
            shutil.rmtree(output_directory) 
        os.mkdir(output_directory) 
        
        dlambda = 0.6
        lambda_obs_min = 3800.0
        lambda_obs_max = 12600.0
        templates_path = "./templates"
        refPath = "/home/aschmitt/data/simu_linemodel/simulm_20160513/simulation_pfswlinemodel_20160513_10spcperbin/refz_bad12.txt"
        refData = loadRef(refPath)
        
        TFDir = "/home/aschmitt/data/simu_linemodel/simulm_20160513/simulation_pfswlinemodel_20160513_10spcperbin"
        catalogPath = "/home/aschmitt/data/simu_linemodel/simulm_20160513/amazed/linecatalogs/linecatalogamazedvacuum_B10E.txt"
        catLines = ctlg.Catalog(catalogPath)

        filenameList = [a[0] for a in refData]
        magList = [a[2] for a in refData]
        tplList = [a[3] for a in refData]
        zList = [a[1] for a in refData]

        print("size of refData = {}".format(len(filenameList))) 
        print("first elt: file = {}, mag = {}, tpl = {}".format(filenameList[0], magList[0], tplList[0]))  
        
        for k in range(len(filenameList)):
            tpl_name = tplList[k]
            z = zList[k]
            
            magRoundedRef = magList[k]
            #loop on the mag uncertainty, yes the mag ref value has been rounded +-0.1 during the creation of the
            magRoundedCheckList = [magRoundedRef-0.05, magRoundedRef-0.025, magRoundedRef, magRoundedRef+0.025, magRoundedRef+0.05]
            filename = filenameList[k]
            #prepare the template spc with the given mag
            tplPath = os.path.join(templates_path, tpl_name)
            tpl = spectrum.Spectrum(tplPath ,'template', snorm=False)

            tfPath = os.path.join(TFDir, "{}_TF.fits".format(filenameList[k]))            
            spcTF = spectrum.Spectrum(tfPath ,'pfs', snorm=False)
            print("tfpath = {}".format(tfPath))
            
            maskOutLines = catLines.getMaskOutsideLines(spcTF.xvect, z)
            
            tplCheckList = []
            lstSqList = []
            for mag in magRoundedCheckList:
                tplTmp = tpl.copy()
                tplTmp.extendWavelengthRangeRed(lambda_obs_max)
                tplTmp.extendWavelengthRangeBlue(lambda_obs_min)
            
                tplTmp.applyLyaExtinction(redshift=z)
                tplTmp.applyRedshift(z)   
#                #Check magFound and mag to be the same
                magFound = tplTmp.setMagIAB(mag)
#                #
#                if np.floor(magFound*10.0)!=np.floor(magRoundedRef*10.0):
#                    print("magRoundedRef = {}, and magFound={}".format(magRoundedRef, magFound))
#                else:
#                    break                
                tplTmp.applyLambdaCrop(lambda_obs_min, lambda_obs_max)
                tplTmp.interpolate(dx=0.1) #high sampling for the synthesis process
                tplTmp.correctZeros()
                tplTmp.interpolate(dx=dlambda)
                
                tplCheckList.append(tplTmp)
                
                lstErr = 0.0
                for kk in range(spcTF.n-1):
                    if maskOutLines[kk]==1.0:
                        lstErr += (spcTF.yvect[kk] - tplTmp.yvect[kk])**2
                print("for k={}, mag={}, lstSq = {}".format(k, mag, lstErr))
                lstSqList.append(lstErr)
            kmin = np.argmin(lstSqList)
            print("found argmin = {}".format(kmin))
            tpl = tplCheckList[kmin]
            
            exported_fits_path = tpl.exportFits(output_directory, "{}_B".format(filename), addNoise=False, exportNoiseSpectrum=False, noiseSigma = 1.5e-19, addNameSuffix = False)
            exported_fits_name = os.path.split(exported_fits_path)[1]
            print("Exported: {}".format(exported_fits_name))
            #stop
        
        print("END of continuum creation script.")
        
    except (KeyboardInterrupt):
        exit()
        
if __name__ == '__main__':
    print "..."
    Main( sys.argv )
    
    
    

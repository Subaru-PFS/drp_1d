# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 10:06:45 2016

@author: aschmitt
"""

import os
import sys
from astropy.io import fits
import argparse
import random
import math
       
import matplotlib.pyplot as plt

import numpy as np
from scipy import interpolate
from scipy.signal import savgol_filter

import spectrum as sp

class PsfTools(object):
    def __init__(self, psftype="instrument"):
        self.psftype = psftype
    
    def applyResolution(self, spc, R=3000.0):
        optionFixed = False
        
        spcModRaw = sp.Spectrum("", stype="empty")
        spcModRaw.setData(spc.xvect, spc.yvect)
        
        #optionnally regrid on a regular grid
        enableRegrid = True
        if enableRegrid:
            spcMod = sp.Spectrum("", stype="empty")
            spcMod.setData(spcModRaw.xvect, spcModRaw.yvect)
            dlRaw = spcModRaw.getResolution()
            
            xvect = np.arange(spcModRaw.xvect[0], spcModRaw.xvect[-1], dlRaw)
            print("Regrid: x first is : {}".format(spcModRaw.xvect[0]))
            print("Regrid: xvect first is : {}".format(xvect[0]))
            print("Regrid: x last is : {}".format(spcModRaw.xvect[-1]))
            print("Regrid: xvect last is : {}".format(xvect[-1]))
            
            spcMod.interpolateOnGrid(xvect)
        else:
            spcMod = spcModRaw
            

        #find input spc sampling rate
        dlambdaSpcMinMax = spcMod.getMaxMinLambdaSteps()
        if np.abs(dlambdaSpcMinMax[0] - dlambdaSpcMinMax[1])<np.abs(dlambdaSpcMinMax[0])*1e-5:
            dlambdaSpc = dlambdaSpcMinMax[0]
            print("INFO: found dlambdaSpc={}".format(dlambdaSpc))
        else:
            print("ERROR: dlambdaSpcMin={}, dlambdaSpcMax={}".format(dlambdaSpcMinMax[0], dlambdaSpcMinMax[1]))
            print("ERROR: input spc doesn't have a fixed klambda grid, aborting...")
            return
        
        if optionFixed:
            #create kernel
            #wlFixedKernel = spcModRaw.xvect[-1]
            wlFixedKernel = 12000.0 #using 12000A fixed wl value
            print("Wl used for fixed kernel = {}".format(wlFixedKernel))
            print("Resolution used for fixed kernel = {}".format(R))
            ker, kerwl = self.createFixedKernel(dlambda=dlambdaSpc, R=R, wl= wlFixedKernel)
            #now apply the convolution kernel        
            yvect_convolved  = np.convolve(spcMod.yvect, ker, mode='same')
            spcMod.yvect = yvect_convolved
            spcMod.name = "{}_psf-r{}".format(spc.name, R)
            print("SPC new resolution name = {}".format(spcMod.name))
        else:
            #loop on spcMod.xvect
            inputFlux = np.array(spcMod.yvect)
            for k, w in enumerate(spcMod.xvect):
                print("Processing wl {}/{}".format(k+1, len(spcMod.xvect)))
                ker, kerwl = self.createFixedKernel(dlambda=dlambdaSpc, R=R, wl= w, enablePlot=False)
                idxZeroKerWl = np.argmin(abs(kerwl-0))
                #print("Info: kernel zero wavelength index = {}".format(idxZeroKerWl))

                sizeKer = len(ker)
                iKerMin = 0
                iKerMax = sizeKer-1
                iFluxMin = k - idxZeroKerWl
                iFluxMax = k + (sizeKer-idxZeroKerWl-1)
                
                #print("Info: iFluxMin = {}, and iFluxMax = {}".format(iFluxMin, iFluxMax))
                #print("Info: iKerMin = {}, and iKerMax = {}".format(iKerMin, iKerMax))
                
                if iFluxMin<0:
                    iKerMin = -iFluxMin
                    iFluxMin = 0
                    
                if iFluxMax>len(inputFlux)-1:
                    iKerMax = sizeKer-( iFluxMax-(len(inputFlux)-1) )-1 
                    iFluxMax = len(inputFlux)-1
                                    
                
                #print("Info: iFluxMin = {}, and iFluxMax = {}".format(iFluxMin, iFluxMax))
                #print("Info: iKerMin = {}, and iKerMax = {}".format(iKerMin, iKerMax))
                spcMod.yvect[k] = np.sum(np.dot(ker[iKerMin:iKerMax], inputFlux[iFluxMin:iFluxMax]))
    
        return spcMod

    def createFixedKernel(self, dlambda, R, wl, enablePlot=True):
        nsigma = 15.0 #total nsigma so that there is nsigma/2 on each side of the gaussian profile
        instrumentResolutionEmpiricalFactor = 230.0/325.0/2.35; 
        sigma = wl/R*instrumentResolutionEmpiricalFactor;        
        npts = int(sigma*nsigma/dlambda)
        print("npts for the kernel is {}".format(npts))
        
        wlvect = np.linspace(0.0, (npts-1)*dlambda, npts) - int(npts/2.0)*dlambda 
        kervect = np.exp(-0.5*wlvect**2/sigma**2)
        sumKervect = np.sum(kervect)
        kervect /= sumKervect
        if enablePlot:
            plt.figure()
            plt.plot(wlvect, kervect)
            plt.title("kernel (fixed: sampling={}, R={}, wavelength={})".format(dlambda, R, wl))
            plt.ylabel("f")
            plt.xlabel("wl (Angstrom)")
            plt.grid()
            
        
        #print("wl support for the kernel is {}".format(wlvect))
        return kervect, wlvect
        
def processSpc(spcpath, stype, enableExport, resolution, dx):
    
    spc = sp.Spectrum(spath = spcpath, stype=stype)
    #spc.applyLambdaCrop(2000, 12000)
    psft = PsfTools(psftype="instrument")
    spcMod = psft.applyResolution(spc, R=resolution)
    spcMod.interpolate(dx=dx)
    spc.plotCompare(spcMod, modellinetype = "-b+")
    
    if enableExport:
        path = os.path.split(spcpath)[0]
        print("export: path = {}".format(path))
        nameWext = os.path.split(spcpath)[1]
        print("export: nameWext = {}".format(nameWext))
        newName = "{}_convolvedR{}".format(os.path.splitext(nameWext)[0], resolution)
        print("export: new name = {}".format(newName))
        if 0:        
            spcMod.exportFits(path, name=newName, addNoise=False, exportNoiseSpectrum=False)
        else:
            fullpath = os.path.join(path, "{}.dat".format(newName))
            spcMod.applyLambdaCrop(1950, 11950)
            spcMod.saveTpl(fullpath)
        
        
def processDir(dirPath, stype, enableExport, resolution, dx):
    fileList = listFiles(dirPath, stype)
    for _file in fileList:
        print("\n")
        print("INFO: processing file: {}".format(_file))
        fullPath = os.path.join(dirPath, _file)
        processSpc(spcpath=fullPath, stype=stype, enableExport=enableExport, resolution=resolution, dx=dx)

def listFiles(f1Directory, stype, exclstring="xrfctgrefardtasrdtrdsatrdtartcomplicatedstring"):
    """
    Scans input directories and build files list whose shape is Nx1
    """
    fileList = []
    if stype=="template":
        ext1 = '.dat'
        ext2 = '.txt'
    else:
        ext1 = '.fits'
        ext2 = '.FITS'
    
 
    for i,file in enumerate(sorted(os.listdir(f1Directory))):
        if (file.endswith(ext1) or file.endswith(ext2)) and exclstring not in file:  
            a1 = str(file)
            #print("file found ={}".format(a1))
            fileList.append(a1)
    return fileList  

def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-i", "--inputPath", dest="inputPath", default="",
                    help="path to the input spectrum or directory to be processed : plotted, etc...")    
    parser.add_argument("-t", "--type", 
                        help="type of spectrum, can be in teh following list \
                        {template, template-w0f2, vvds, pfs, pfs2, muse, hplana, euclid, pfs2reech, euclid_sim_noise, euclid_sim_siroct2016_flux, euclid_sim_siroct2016_noise}",  
                        dest="spcType",     
                        default="") 
    parser.add_argument("-e", "--export", dest="export", action='store_true',
                    help="enable export of the convolved spectrum")
                    
                    
    parser.add_argument("-r", "--resolution", dest="resolution", default="5000",
                    help="instrument resolution")
    parser.add_argument("-d", "--deltalambda", dest="deltalambda", default="0.6",
                    help="wavelength sampling")
                    
    options = parser.parse_args()
                    
    resolution = float(options.resolution)
    dx = float(options.deltalambda)
    print("INFO: using r={}".format(resolution))
    print("INFO: using dlambda={}".format(dx))

    if os.path.exists(options.inputPath):
        if os.path.isdir(options.inputPath):
            processDir(dirPath=options.inputPath, stype=options.spcType, enableExport=options.export,
                   resolution=resolution, dx=dx)
        elif os.path.isfile(options.inputPath):
            processSpc(spcpath=options.inputPath, stype=options.spcType, enableExport=options.export,
                   resolution=resolution, dx=dx)
        else:
            print("ERROR: input was not a file, nor a directory... aborting")
    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print("psftools")
    Main( sys.argv )
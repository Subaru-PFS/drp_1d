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
        

def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-s", "--spc", dest="spcPath", default="",
                    help="path to the fits spectrum to be plotted")
    parser.add_argument("-e", "--export", dest="export", action='store_true',
                    help="enable export of the convolved spectrum")
    options = parser.parse_args()
                    

    if os.path.exists(options.spcPath):
        spc = sp.Spectrum(spath = options.spcPath)
        psft = PsfTools(psftype="instrument")
        spcMod = psft.applyResolution(spc, R=300.0)
        spcMod.interpolate(dx=0.6)
        spc.plotCompare(spcMod, modellinetype = "-b+")
        
        if options.export:
            path = os.path.split(options.spcPath)[0]
            print("export: path = {}".format(path))
            nameWext = os.path.split(options.spcPath)[1]
            print("export: nameWext = {}".format(nameWext))
            newName = "{}_convolved".format(os.path.splitext(nameWext)[0])
            print("export: new name = {}".format(newName))
            spcMod.exportFits(path, name=newName, addNoise=False, exportNoiseSpectrum=False)

    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print("Spectrum")
    Main( sys.argv )
# -*- coding: utf-8 -*-
#/usr/bin/env python
"""
Created on Sat Jul 25 11:44:49 2015

@author: aschmitt
"""
import os
import sys
from astropy.io import fits
import pyfits
import optparse
import random
import math

from bokeh.plotting import figure, output_file, show
        
import matplotlib.pyplot as pp
import numpy as np
from scipy import interpolate

   
class Spectrum(object):
    def __init__(self, spath, stype='undefspc', snorm=False, label=""):
        self.logTagStr = "Spectrum"
        self.spath = spath
        self.name = os.path.basename(spath)
        self.label = label
        self.stype = stype
        if self.stype == "" or self.stype == "-1":
            self.stype = 'undefspc'
        print("type = {0}".format(stype))
        print("norm = {0}".format(snorm))
        self.snorm = snorm
        self.n = -1
        self.xvect = -1
        self.yvect = -1
        self.ysum = 0
        self.forcePlotXIndex = False
        if(self.stype == 'template'):
            self.loadTpl() 
        elif(self.stype == 'pfs2reech'):
            self.loadpfs2reech() 
        elif(self.stype == 'pfs2'):
            self.loadpfs2() 
        elif(self.stype == 'hplana'):
            self.loadpfs2() 
        elif(self.stype == 'vvds'):
            self.loadvvds() 
        elif(self.stype == 'pfs'):
            self.loadpfs() 
        elif(self.stype == 'muse'):
            self.loadmuse() 
        else:
            self.load()
               
        if self.snorm:
            print("WARNING: snorm = {0}".format(snorm))
            for x in range(0,self.n):
                self.yvect[x] /= self.ysum/self.n
            self.ysum = 0.0
            for x in range(0,self.n):
                self.ysum += self.yvect[x]

    def copy(self):
        scopy = Spectrum(self.spath, self.stype, self.snorm)
        return scopy
        
    def load(self):
        hdulist = fits.open(self.spath) 
        print("\nHDULIST = \n{}\n".format(hdulist))
        try:
            sciheader = hdulist[1].header
        except: 
            sciheader = hdulist[0].header
            
        #print("header = {0}".format(sciheader))
        n = sciheader["NAXIS"]
        #n1 = sciheader["NAXIS1"]
        
        try:
            n2 = sciheader["NAXIS2"]
        except:
            print("INFO: USING MUSE spectrum format...")
            self.stype='muse' #euclid and muse splitted
            self.loadmuse()
            return
                    
        #print("n1 = {0}".format(n1))
        #print("n2 = {0}".format(n2))
        if self.stype == 'undefspc' and n2 > 2:
            self.stype='pfs'
            self.loadpfs()
        elif self.stype == 'undefspc' and n2 < 2 and n<2:
            self.stype='pfs2reech'
            self.loadpfs2reech()
        elif self.stype == 'undefspc' and n2 < 2:
            self.stype='vvds'
            self.loadvvds()

    def loadpfs(self):
        hdulist = fits.open(self.spath) 

        try:
            scidata = hdulist[1].data
        except:
            scidata = hdulist[0].data
            
        self.n = scidata.shape[0]
        if self.n<3:
            self.n = scidata.shape[1]
        #print('{0} - n = {1}'.format(self.logTagStr, self.n))
    
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            #print scidata[x]
            self.xvect[x] = scidata[x][0]
            self.yvect[x] = scidata[x][1]
            self.ysum += self.yvect[x]
    
    #same as pfs but with xaxis in nm instead of Angstrom
    def loadpfs2(self):
        hdulist = fits.open(self.spath) 

        try:
            scidata = hdulist[1].data
        except:
            scidata = hdulist[0].data
            
        self.n = scidata.shape[0]
        if self.n<3:
            self.n = scidata.shape[1]
        #print('{0} - n = {1}'.format(self.logTagStr, self.n))
    
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            self.xvect[x] = scidata[x][0]*10.0
            self.yvect[x] = scidata[x][1]
            self.ysum += self.yvect[x]
            
    def loadpfs2reech(self):
        """
        this type corresponds to pfs spectra reech. by Sara
        """
        hdulist = fits.open(self.spath) 

        try:
            scidata = hdulist[1].data
        except:
            scidata = hdulist[0].data
            
        self.n = scidata.shape[0]
        if self.n<3:
            self.n = scidata.shape[1]
        #print('{0} - n = {1}'.format(self.logTagStr, self.n))

        scidata = hdulist[0].data
        sciheader = hdulist[0].header
        
        if 'LAMBDA' in sciheader["CTYPE1"]:
            xaxisStart = sciheader["CRVAL1"]
            xaxisRes = sciheader["CDELT1"]    
    
        #print sciheader
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            self.xvect[x] = (xaxisStart + (x)*xaxisRes) #warning: nm unit for this case
            #print('scidata[0]= {}'.format(scidata[0]))
            self.yvect[x] = scidata[x]
            self.ysum += self.yvect[x]
        
    def loadvvds(self):
        hdulist = fits.open(self.spath) 

        scidata = hdulist[0].data
        sciheader = hdulist[0].header
        
        if 'LAMBDA' in sciheader["CTYPE1"]:
            xaxisStart = sciheader["CRVAL1"]
            xaxisRes = sciheader["CDELT1"]
            xaxisRefPix = sciheader["CRPIX1"]
            
            
        self.n = scidata.shape[0]
        if self.n<3:
            self.n = scidata.shape[1]
        #print('{0} - n = {1}'.format(self.logTagStr, self.n))
    
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            self.xvect[x] = xaxisStart + (x+1-xaxisRefPix)*xaxisRes
            self.yvect[x] = scidata[0][x]
            self.ysum += self.yvect[x]


    def loadmuse(self):
        hdulist = fits.open(self.spath) 

        scidata = hdulist[0].data
        sciheader = hdulist[0].header
        
        xaxisStart = sciheader["CRVAL1"]
        xaxisRes = sciheader["CDELT1"]
        xaxisRefPix = sciheader["CRPIX1"]
            
            
        self.n = scidata.shape[0]
        if self.n<3:
            self.n = scidata.shape[1]
        #print('{0} - n = {1}'.format(self.logTagStr, self.n))
    
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            self.xvect[x] = xaxisStart + (x+1-xaxisRefPix)*xaxisRes
            self.yvect[x] = scidata[x]
            self.ysum += self.yvect[x]            
            
    def loadhplana(self):
        print("Loading hplana fits")
        hdulist = fits.open(self.spath) 
        print hdulist
        
        scidata = hdulist[0].data
        sciheader = hdulist[0].header
        
        xaxisStart = sciheader["CRVAL1"]
        print("xaxisstart = {}".format(xaxisStart))
        xaxisRes = sciheader["CDELT1"]
        print("xaxisRes = {}".format(xaxisRes))
        xaxisRefPix = sciheader["CRPIX1"]
        print("xaxisRefPix = {}".format(xaxisRefPix))
            
            
        self.n = scidata.shape[1]
        #print('{0} - n = {1}'.format(self.logTagStr, self.n))
    
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            self.xvect[x] = xaxisStart + (x+1-xaxisRefPix)*xaxisRes
            self.yvect[x] = scidata[0][x]
            self.ysum += self.yvect[x]
            
                
    def loadmuseraw(self):
        print("Loading muse (raw, non splitted) fits")
        hdulist = fits.open(self.spath) 
        print hdulist

        if 1: #loading data  
            print("Loading DATA")
            scidata = hdulist[1].data
            sciheader = hdulist[1].header
        else: #loading noise        
            print("Loading NOISE")
            scidata = hdulist[2].data
            sciheader = hdulist[2].header
            
        xaxisStart = sciheader["CRVAL1"]
        print("xaxisstart = {}".format(xaxisStart))
        xaxisRes = sciheader["CDELT1"]
        print("xaxisRes = {}".format(xaxisRes))
        xaxisRefPix = sciheader["CRPIX1"]
        print("xaxisRefPix = {}".format(xaxisRefPix))
            
            
        self.n = scidata.shape[0]
        print('{0} - n = {1}'.format(self.logTagStr, self.n))
    
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            self.xvect[x] = xaxisStart + (x+1-xaxisRefPix)*xaxisRes
            self.yvect[x] = scidata[x]
            self.ysum += self.yvect[x]
    
    def loadTpl(self):
        filename = self.spath
        wave = []
        flux = []
        f = open(filename)
        for line in f:
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                #print lineStr
                data = lineStr.split(" ")
                data = [r for r in data if r != '']
                #print len(data)
                if(len(data) >=2):
                    # fill the list                
                    #print data[0]
                    #print data[1]
                    wave.append(float(data[0]))
                    flux.append(float(data[1]))
                else:
                    #try to separate the values with \t instead...
                    data = lineStr.split("\t")
                    data = [r for r in data if r != '']
                    #print len(data)
                    if(len(data) >=2):
                        # fill the list                
                        #print data[0]
                        #print data[1]
                        wave.append(float(data[0]))
                        flux.append(float(data[1]))
        f.close()
        self.n = len(wave)
        #print('len wave = {0}'.format(self.n))
        #---- default xaxis index array
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0    
        for x in range(0,self.n):
            self.xvect[x] = wave[x]
            self.yvect[x] = flux[x]
            self.ysum += self.yvect[x]   
            
    def saveTpl(self, outputfullpath):
        filename = outputfullpath
        f = open(filename, "w")

        f.write("# lambda\tflux\n")
        for x in range(0,self.n):        
            f.write(" " + str(self.xvect[x]) + " " + str(self.yvect[x]) + "\n")
        f.close()
        
        
    def __str__(self):
        a = "\nSpectrum: {0}\n".format(self.name)
        a = a + ("    type = {0}\n".format(self.stype))
        a = a + ("    n = {0}\n".format(self.n))
        a = a + ("    dlambda = {0}\n".format(self.getResolution()))
        a = a + ("    lambda min = {}, lambda max = {}\n".format(self.getWavelengthMin(), self.getWavelengthMax()))
        a = a + ("    flux min = {}, flux max = {}\n".format(self.getFluxMin(),self.getFluxMax()))
        a = a + ("    flux std = {}, flux std 6000_8000 = {}\n".format(np.std(self.yvect), self.GetFluxStd6000_8000()))
        a = a + ("    magI = {}\n".format(self.getMagIAB()))
        
        a = a + ("\n")
        
        return a
        
    def plot(self, saveFullDirPath=""):
        #pp.plot(self.xvect, self.yvect, "x-")
        pp.plot(self.xvect, self.yvect)

        pp.grid(True) # Affiche la grille
        #pp.legend(('cos','sin'), 'upper right', shadow = True)
        if not self.forcePlotXIndex:
             pp.xlabel('angstrom')
        else:
            pp.xlabel('index')
        pp.ylabel('y')
        pp.title(self.name) # Titre
        if not saveFullDirPath == "":
            saveFullPath = os.path.join(saveFullDirPath, self.name + str(".png"))
            pp.savefig(saveFullPath) # sauvegarde du fichier ExempleTrace.png
        pp.show()
        print '\n'
        
    def plotCompare(self, other_spc, amplitude = 1.0, modellinetype = "-bo", exportPath=""):
        fig = pp.figure( "spectrumview", figsize=(15,11))
        
        lbl = self.label
        if lbl=="":
            lbl = self.name
        spcColor = '0.6'
        pp.plot(self.xvect, self.yvect,  color=spcColor, label=self.name)
        yother = [a*amplitude for a in other_spc.yvect]
        lbl = other_spc.label
        if lbl=="":
            lbl = other_spc.name
        pp.plot(other_spc.xvect, yother, modellinetype, label=other_spc.name)

        pp.grid(True) # Affiche la grille
        pp.legend()
        if not self.forcePlotXIndex:
             pp.xlabel('Angstrom')
        else:
            pp.xlabel('index')
        pp.ylabel('Flux')
        pp.title(self.name) # Titre
        if exportPath=="":
            pp.show()
        else:
            pp.savefig(exportPath)
            pp.close()
            pp.clf()
            
        print '\n'
        
    def exportPlotCompareBokeh(self, other_spc, amplitude = 1.0, fname = "spectrum_comparison"):
        # output to static HTML file
        output_file(fname + ".html", title=fname)

        ylabel = 'flux'
        if not self.forcePlotXIndex:
            xlabel = 'angstrom'
        else:
            xlabel = 'index'
            
        # create a new plot with a title and axis labels
        TOOLS="resize,crosshair,pan,wheel_zoom,box_zoom,reset,box_select,lasso_select"
        TOOLS="resize,crosshair,pan,wheel_zoom,box_zoom,reset"
        p = figure(tools=TOOLS, title=fname, x_axis_label= xlabel, y_axis_label=ylabel)

        # add a line renderer with legend and line thickness
        p.line(self.xvect, self.yvect, legend="Spectrum", line_width=1, line_color="black")
        yother = [a*amplitude for a in other_spc.yvect]
        p.line(other_spc.xvect, yother, legend="Model", line_width=2, line_color="blue", line_dash=[4, 4])

        # show the results
        show(p)
    
    def getResolution(self):
        dl = 0.0;
        for x in range(0,self.n-1):
            dl += self.xvect[x+1]-self.xvect[x]
        dl /= self.n
        return dl
        
    def getShiftedX(self, z):
        coeff= 1+z
        xvect = range(0,self.n)
        for x in range(0,self.n):
            xvect[x] = self.xvect[x]*coeff
        return xvect
    
    def getWeightedY(self, a):
        yvect = range(0,self.n)
        for x in range(0,self.n):
            yvect[x] = self.yvect[x]*a
        return yvect
        
    def getWeightedFlux(self, xval, a):
        if xval < self.xvect[0]:
            return self.xvect[0]*a
        elif xval > self.xvect[self.n-1]:
            return self.xvect[self.n-1]*a
        else:
            for x in range(0,self.n):
                if  xval > self.xvect[x]:
                    return self.yvect[x]*a
                    
    def getFluxMin(self):
        return min(self.yvect)                    
    def getFluxMax(self):
        return max(self.yvect)
        
    def GetFluxStd6000_8000(self):
        imin = self.getWavelengthIndex(6000)
        imax = self.getWavelengthIndex(8000)
        if imin == 0 or imax == 0:
            return -1.0
        return np.std(self.yvect[imin:imax])
    
    def getWavelengthMin(self):
        return self.xvect[0]
    def getWavelengthMax(self):
        return self.xvect[len(self.xvect)-1]    
    def getWavelengthIndex(self, wl):
        i=0
        for x in range(0,self.n-1):
            if wl >= self.xvect[x]:
                i=x
            if wl < self.xvect[x+1]:
                break
        return i
        
        
    def smoothGaussian(self,signal,degree=5):  
        window=degree*2-1  
        weight=np.array([1.0]*window)  
        weightGauss=[]  
        for i in range(window):  
            i=i-degree+1  
            frac=i/float(window)  
            gauss=1/(np.exp((4*(frac))**2))  
            weightGauss.append(gauss) 
        weight=np.array(weightGauss)*weight  
        smoothed=[0.0]*(len(signal)-window)  
        for i in range(len(smoothed)):  
            smoothed[i]=sum(np.array(signal[i:i+window])*weight)/sum(weight) 
        
        dec = len(signal)-len(smoothed)
        decInit = int(dec/2.0)
        decEnd = dec-decInit
        #print dec
        smoothed2=[0.0]*(len(signal)) 
        for i in range(0,decInit):
            smoothed2[i]=smoothed[0]
        for i in range(len(smoothed)):  
            smoothed2[i+decInit]=smoothed[i]
        for i in range(len(smoothed), len(signal)-decInit):  
            smoothed2[i+decInit]=smoothed[len(smoothed)-1]
        return smoothed2 
        
    def extendWavelengthRangeRed(self, wavelengthSup):
        i1 = self.n-2;
        i2 = self.n-1;
        x1 = self.xvect[i1];
        y1 = self.yvect[i1];
        x2 = self.xvect[i2];
        y2 = self.yvect[i2];
        xstep = x2-x1
        ystep = y2-y1  
        #use a poly1 fit
        NRangeforfit = 200
        xv = self.xvect[self.n-NRangeforfit:self.n-1]
        yv = self.yvect[self.n-NRangeforfit:self.n-1]       
        coefficients1 = np.polyfit(xv,yv,1)
        print("linfit coeff.: 0 = {0}, 1 = {1}".format(coefficients1[0], coefficients1[1]))
        moyenne = np.mean(yv)
        
        
        k = 1
        x = x2
        while(x<wavelengthSup and k<1e6):
            x = x2 + k*xstep
            y = coefficients1[1] + coefficients1[0]*x #y2# + k*ystep
            #y = moyenne
            print("INFO (extension Red): new sample at x = {} and y = {}".format(x, y))            
            self.xvect.append(x)
            self.yvect.append(y)
            self.n = self.n + 1
            k=k+1
        if k>=1e6:
            print("WARNING: stopping extension, max iteration reached")

    def extendWavelengthRangeBlue(self, wavelengthInf):
        i1 = 0;
        i2 = 1;
        x1 = self.xvect[i1];
        y1 = self.yvect[i1];
        x2 = self.xvect[i2];
        y2 = self.yvect[i2];
        xstep = x2-x1
        ystep = y2-y1        
        #use a poly1 fit
        NRangeforfit = 200
        xv = self.xvect[0:NRangeforfit]
        yv = self.yvect[0:NRangeforfit]       
        coefficients1 = np.polyfit(xv,yv,1)
        print("linfit coeff.: 0 = {0}, 1 = {1}".format(coefficients1[0], coefficients1[1]))
        moyenne = np.mean(yv)
        
        k = 1
        x = x1
        while(x>wavelengthInf and k<1e6):
            x = x1 - k*xstep
            y = coefficients1[1] + coefficients1[0]*x#y1 #-k*ystep
            #y = moyenne
            
            # custom profile
            if x<912:
                y=0
#            else:
#                x1=2238.0
#                y1=9.11*1e-19
#                x0=912.0
#                y0=0.0
#                a = (y1-y0)/(x1-x0)
#                b = -a*x0
#                y = b + a*x#y1 #-k*ystep
            
            print("INFO (extension Blue: new sample at x = {} and y = {}".format(x, y))
            xvect_tmp = []
            yvect_tmp = []
            for i in range(self.n):
                xvect_tmp.append(self.xvect[i])
                yvect_tmp.append(self.yvect[i])
            self.xvect.append(xvect_tmp[self.n-1])
            self.yvect.append(yvect_tmp[self.n-1])
            self.xvect[0] = x
            self.yvect[0] = y
            
            for i in range(self.n):
                self.xvect[i+1]=xvect_tmp[i]
                self.yvect[i+1]=yvect_tmp[i]
            self.n = self.n + 1
            k=k+1
        if k>=1e6:
            print("WARNING: stopping extension, max iteration reached")
            
    def interpolate(self, dx=1.0):
        x = np.copy(self.xvect)
        y = np.copy(self.yvect)
        f = interpolate.interp1d(x, y)

        self.xvect = np.arange(x[0], x[self.n-1], dx)
        self.yvect = f(self.xvect) 
        self.n = len(self.yvect)
        
    def applyWeight(self, w):
        for x in range(0,self.n):
            self.yvect[x] = self.yvect[x]*w
    
    def applyRedshift(self, z):
        coeff= 1+z
        for x in range(0,self.n):
            self.xvect[x] = self.xvect[x]*coeff
    
    def applyLambdaCrop(self, lambdaMin, lambdaMax):
        imin = 0
        imax = self.n-1
        for k in range(1,self.n):
            if self.xvect[k-1]<=lambdaMin and self.xvect[k]>=lambdaMin:
                imin = k
            if self.xvect[k-1]<=lambdaMax and self.xvect[k]>=lambdaMax:
                imax = k
        old_xvect = self.xvect
        old_yvect = self.yvect
            
        print("imin = {}, imax = {}".format(imin, imax))
        self.n = imax-imin+1
        self.xvect = range(0,self.n)
        self.yvect = range(0,self.n)
        self.ysum = 0.0
        for x in range(0,self.n):
            self.xvect[x] = old_xvect[x+imin]
            self.yvect[x] = old_yvect[x+imin]
            self.ysum += self.yvect[x]
    
            
    def setMagIAB(self, magIAB):
        dMThreshold = 0.1
        maxIterations = 10000
        mag = self.getMagIAB()     
        
        it = 0
        coeffStep = 1.+dMThreshold*10.0
        while abs(mag-magIAB) > dMThreshold and it<maxIterations:
            if mag>magIAB:
                self.applyWeight(coeffStep)
            else:
                self.applyWeight(1/coeffStep)
            mag = self.getMagIAB() 
            #print("current Mag = {}".format(mag))
            it +=1
            if it>20:
                coeffStep = 1.+dMThreshold
        print("Mag found (n={} iterations) is {}".format(it, mag))


        
    def getMagIAB(self):
        lambda_min=6000
        lambda_max=9800
        imin = self.getWavelengthIndex(lambda_min)
        imax = self.getWavelengthIndex(lambda_max)

        
        y = self.getFNU()[imin:imax]
        #y = self.yvect[imin]*self.xvect[imin]*self.xvect[imin]
        #f = np.mean(y)/1e-29/3e18
        f = np.mean(y)
        
        #sumF = 8.31*1e-9#2.55*10e-20 #Vega, band I
        #print("meanF = {}".format(meanF))
        #refF = 3631.0#2550*1e-23#8.31*1e-9#2.55*10e-20 #Vega, band I
        #refM = 48.6
        mag = -5.0/2.0*np.log10(f)-48.6
        #Check: V. Lebrun's numbers F = 1e-18 at mag 23.5
        return mag
    
    def getFNU(self):
        c_cm_s = 3e18
        yvect = np.array(self.yvect)
        
        for x in range(self.n):
            #self.yvect[x] = self.yvect[x]*3.34*1e4*self.xvect[x]**2
            yvect[x] = self.yvect[x]/c_cm_s*(self.xvect[x]**2)
            
        return yvect

    def exportFits(self, path="", name="", addNoise=False, exportNoiseSpectrum=False, noiseSigma = 1.5e-19):
        """
        write spectrum into a new fits file (useful for model csv input spectra)
        """
        if name == "":
            name = "spectrum_exported"
        noisename = "{}_ErrF".format(name)
        if addNoise:
            name = "{}_F".format(name)
        else:
            name = "{}_TF".format(name)
        destfileout = os.path.join(path, "{}.fits".format(name))
        print("Spectrum exporting to fits: {}".format(destfileout))
        if os.path.exists(destfileout):
            print("Spectrum deleting existing fits: {}".format(destfileout))
            os.remove(destfileout)

        anoise = np.ones((len(self.yvect)))
        snr_meanFlux = np.mean(self.yvect)
        snr_rmsNoise = 0.0
        if addNoise or exportNoiseSpectrum:
            #noise
            #knoise = 0.1 * np.std(self.yvect)
            sigma_from_pfs_simu = 1.5e-19 #corresponds to 3h exposure time
            #noiseSigma = sigma_from_pfs_simu 
            sigma_from_vvds_typical = 1.5e-19 #corresponds approx. to SNR=7 at Mag=23.5
            print("using Noise STD fixed value = {}".format(noiseSigma))            
            
            nvalues = len(self.yvect)
            for k in range(nvalues):
                mu = 0.0
                sigma = np.sqrt(noiseSigma**2)#+self.yvect[k]**2) todo: add the poisson component
                anoise[k] =  sigma
                if addNoise:
                    #print("Adding Noise to spectrum".format()) 
                    noise = random.gauss(mu, sigma)
                    snr_rmsNoise += noise*noise
                    self.yvect[k] = self.yvect[k] + noise
            if addNoise:
                snr_rmsNoise = np.sqrt(snr_rmsNoise/(nvalues-1))
                snr = snr_meanFlux/snr_rmsNoise
                print("SNR for the exported spectrum is: {}".format(snr))
        
                  
        a1 = self.xvect
        #print("col1={}".format(a1))
        a2 = self.yvect
        #print("col2={}".format(a2))
        col1 = pyfits.Column(name='wave', format='E', array=a1)
        col2 = pyfits.Column(name='flux', format='E', array=a2)
        cols_spc = pyfits.ColDefs([col1, col2])
 
        #INFO: This shortcut will automatically create a minimal primary HDU with no data and prepend it to the table HDU to create a valid FITS file.
        # from https://pythonhosted.org/pyfits/ 
        hdu_new = pyfits.BinTableHDU.from_columns(cols_spc)
        hdu_new.writeto(destfileout)
        
        if exportNoiseSpectrum:
            
            destfileoutnoise = os.path.join(path, "{}.fits".format(noisename))
            print("Noise exporting to fits: {}".format(destfileoutnoise))
            if os.path.exists(destfileoutnoise):
                print("Spectrum deleting existing fits: {}".format(destfileoutnoise))
                os.remove(destfileoutnoise)   
            
            colnoise = pyfits.Column(name='noise', format='E', array=anoise)
            cols_noise = pyfits.ColDefs([col1, colnoise])
            hdu_noise = pyfits.BinTableHDU.from_columns(cols_noise)
            hdu_noise.writeto(destfileoutnoise)
            
        print ""
    
        return destfileout
            

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./spectum.py -s """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-s", u"--spc", help="path to the fits spectrum to be plotted",  dest="spcPath", default="")
    parser.add_option(u"-t", u"--type", help="type of spectrum, can be in teh following list {template, vvds, pfs, pfs2, muse, hplana, euclid, pfs2reech}",  dest="spcType", default="")
    parser.add_option(u"-e", u"--export", help="export to fits format",  dest="export", default="no")

    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        print('using full path: {0}'.format(options.spcPath))
        s = Spectrum(options.spcPath, options.spcType)
        #s.applyRedshift(0.25)
        
        
        if options.export == "yes":
            #s.applyLambdaCrop(7500, 9000)
            #s.applyWeight(1e-18)
            #s.setMagIAB(20)
            path = os.path.split(options.spcPath)[0]
            nameWext = os.path.split(options.spcPath)[1]
            s.exportFits(path, name=os.path.splitext(nameWext)[0], addNoise=False, exportNoiseSpectrum=False)
        print(s) 
        s.plot()
    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print "Spectrum plotting"
    if 1:
        Main( sys.argv )
           
    if 0: #deprecated, but not fully covered by the command line args functionnality...
        print "Spectrum plotting"
        # plot single spectrum
        if 0:
            path = "/home/aschmitt/data/pfs/pfs_lbg/lbgabs_1K_2z3_20J22.5"
            #name = "EZ_fits-W-ErrF_9.fits"
            #name = "EZ_fits-W-F_9.fits"
            name = "EZ_fits-W-TF_445.fits"
            
            path = "/home/aschmitt/data/pfs/pfs_reallyjustline/reallyjustlinecont_1k_0.5z1.8_0.1A100/"
            #name = "EZ_fits-W-TF_929.fits"
            
            path = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/../1D"
            name = "sc_020086397_F02P016_vmM1_red_31_1_atm_clean.fits"
              
            #path = "/home/aschmitt/data/pfs/all_fromSara/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100"
            #path = "/home/aschmitt/data/pfs/all_fromSara/reech_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100"        
            #name = "SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version0.fits"
            
            #path = "/home/aschmitt/data/pfs/reech_fromSara/PFS_data_FITS_reechantillonne_noRMSE_trueflux_FILES/trueFlux/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100/"
            #name = "SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version929.fits"
    
            path = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/fits"
            name = "NEW_Im_extended_blue.fits"
            
            path = "/home/aschmitt/data/pfs/pfs_testsimu_20151105/470026900000130.24-0.426_20.3_20.4vacLine"
            name = "EZ_fits-W-TF_0.fits"
            
            path = "/home/aschmitt/data/hplana/Residual_fits"
            name = "0278-51900-0151.fits"
    
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s = Spectrum(spath, "hplana")
            print(s) 
            s.plot()
        
        # plot single template
        if 0:
            path = "/home/aschmitt/code/python/simulation_perf_matrix_linemodel/templates"
            name = "NEW_Im_extended.dat"
            #path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/galaxy"
            #name = "sadataExtensionData.dat"
            
            #path = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/output/sc_020103374_F02P021_M1_red_49_1_atm_clean"
            #name = "linemodelsolve.linemodel_spc_extrema_0.csv"
    
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s = Spectrum(spath, 'template', snorm=True)
            print(s)
            s.interpolate()
            #s.exportFits()
            s.plot()   
            
        # plot single template, interp. normalized, exported as fits
        if 0:
            path = "/home/aschmitt/code/python/simulation_perf_matrix_linemodel/templates"
            name = "spectrum_exported.fits"
            #path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/galaxy"
            #name = "sadataExtensionData.dat"
            
            #path = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/output/sc_020103374_F02P021_M1_red_49_1_atm_clean"
            #name = "linemodelsolve.linemodel_spc_extrema_0.csv"
    
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s = Spectrum(spath)
            print(s)
            s.plot()
        
        # plot model
        if 0:
            path = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/output/sc_020103374_F02P021_M1_red_49_1_atm_clean"
            name = "linemodelsolve.linemodel_spc_extrema_0.csv"
    
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s = Spectrum(spath, 'template')
            print(s) 
            s.plot()
            
                
        # extend and plot single template
        if 0:
            path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/emission"
            name = "NEW_Im_extended_blue.dat"
            #path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/galaxy"
            #name = "sadataExtensionData.dat"
            
            path = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/output/sc_020103374_F02P021_M1_red_49_1_atm_clean"
            name = "linemodelsolve.linemodel_spc_extrema_0.csv"
            
            path = "/home/aschmitt/data/pfs/pfs2_simu20151118_jenny/amazed/templates/ExtendedGalaxyEL2/galaxy"
            name = "EW_SB2extended.dat"
    
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s = Spectrum(spath, 'template')
            s.extendWavelengthRangeRed(13000)
            s.extendWavelengthRangeBlue(100)
            outputpath = "/home/aschmitt/data/pfs/pfs2_simu20151118_jenny/amazed/templates/MoreExtendedGalaxyEL2/galaxy"
            soutputpath = os.path.join(outputpath,name)
            s.saveTpl(soutputpath)
            print(s) 
            s.plot()
            
        # compare templates
        if 0:          
            path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/emission"
            name = "NEW_Im_extended_blue.dat"
            #name = "NEW_Im_extended.dat" #1
            #name = "NEW_Sbc_extended.dat" #1        
            #name = "Scd.txt" #1
            #name = "StarBurst1.txt" #2
            #name = "StarBurst2.txt" #2
            #name = "StarBurst3.txt" #2
            
            path = "/home/aschmitt/data/pfs/pfs2_simu20151118_jenny/amazed/templates/ExtendedGalaxyEL3/emission"
            #name = "zcosmos_red.txt" #1
            #name = "s0dataExtensionData.dat"    #3    
            #name = "BulgedataExtensionData.dat"  #3      
            name = "EllipticaldataExtensionData.dat" #3
            
            name = "EW_SB2extended.dat" #2
    
            
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s1 = Spectrum(spath, 'template', snorm=True)
            print(s1)
            #s.plot()
        
            #path = "/home/aschmitt/gitlab/amazed/bin"
            #name = "template_fine.txt"
            path = "/home/aschmitt/data/pfs/pfs2_simu20151118_jenny/amazed/templates/ExtendedGalaxyEL3/galaxy" 
            name = "NEW_E_extendeddataExtensionData.dat"
            #name = "BulgedataExtensionData.dat"  #3
            #name = "sadataExtensionData.dat"
            #name = "EdataExtensionData.dat"
    
            #path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/emission"
            #name = "Scd.txt" #1
            
            #name = "StarBurst2.txt" #2
            
            #path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/galaxy"
            #name = "BulgedataExtensionData.dat"
            #path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/galaxy"
            #name = "EllipticaldataExtensionData.dat"
            
            #path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/galaxy"
            #name = "EW_SB2extended.dat"
            
            
            #path = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/emission"
            #name = "NEW_Im_extended_blue_continuum.txt"
    
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s2 = Spectrum(spath, 'template', snorm=True)
            #s2.yvect =  s2.smoothGaussian(s2.yvect,degree=10)
            print(s2)
            
            s1.plotCompare(s2, 1.0)
            #s1.plotCompare(s2, s1.ysum/s2.ysum)
            
        # compare spectrum
        if 0:
            path = "/home/aschmitt/data/pfs/pfs_reallyjustline/reallyjustlinecont_1k_0.5z1.8_0.1A100/"
            name = "EZ_fits-W-F_622.fits"
            path = "/home/aschmitt/data/pfs/pfs_reallyjustline/reallyjustlinecont_1k_0.5z1.8_0.1A100/"
            name = "EZ_fits-W-F_0.fits"
            path = "/home/aschmitt/data/pfs/all_fromSara/reech_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100"        
            name = "SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version0.fits"
            path = "/home/aschmitt/data/pfs/pfs_reallyjustline/reallyjustlinecont_1k_0.5z1.8_0.1A100/"
            name = "FILTERED_REECH_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version0.fits"
    
            path = "/home/aschmitt/gitlab/amazed/bin/"
            name = "spectrum.fits"     
            #name = "spectrum4linefit.fits"     
         
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s1 = Spectrum(spath)
            print(s1) 
            #s1.plot() 
            
            #path = "/home/aschmitt/data/pfs/reech_fromSara/PFS_data_FITS_reechantillonne_noRMSE_trueflux_FILES/trueFlux/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100/"
            #name = "SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version929.fits"
            path = "/home/aschmitt/data/pfs/all_fromSara/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100/SET_reallyjustlinecont_1k_0.5z1.8_0.1A100"
            name = "SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version0.fits"
            
            path = "/home/aschmitt/data/pfs/pfs_reallyjustline/reallyjustlinecont_1k_0.5z1.8_0.1A100/"
            name = "FILTERED_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version0.fits"
            
            path = "/home/aschmitt/gitlab/amazed/bin/"
            name = "model.fits"
            #name = "model4linefit.fits"
            spath = os.path.join(path,name)
            print('using full path: {0}'.format(spath))
            s2 = Spectrum(spath)
            print(s2) 
            
            s1.plotCompare(s2, 1.0, modellinetype = "bo-")    
               
            
# -*- coding: utf-8 -*-
"""
@author: aschmitt
"""

import os
import sys
import re
import random

import matplotlib as mpl
mpl.use('Qt5Agg')

import matplotlib.pyplot as pp
import numpy as np

import argparse

class ModelResult(object):
    def __init__(self, modelpath):
        self.logTagStr = "Model"
        
        self.name = os.path.basename(modelpath)
        
        self.modelpath = modelpath
        self.n = 0
        
        self.load()
        
    def load(self):
        filename = self.modelpath
        
        name = []
        linelambda = []
        force = []
        linetype = []
        
        lineid = []
        amplitude = []
        error = []
        errorfit = []
                
        f = open(filename)
        for line in f:            
            #print("line = {}".format(line))
            lineStr = line.strip()
            if not lineStr.startswith('#'):
                data = re.split("\t| ", lineStr)
                data = [r for r in data if r != '']
                #print("len(data) = {}".format(len(data)))
                if(len(data) >=8):  
                    name.append(data[2])
                    linelambda.append(float(data[4]))
                    force.append(data[1])
                    linetype.append(data[0])
                    
                    lineid.append(float(data[3]))
                    amplitude.append(float(data[5])) 
                    error.append(float(data[6])) 
                    errorfit.append(float(data[7]))
            if lineStr.startswith('#lya asymfit params = {'):
                data = re.split("widthCoeff=", lineStr)
                data = re.split(", alpha=", data[1])
                self.lya_widthCoeff = float(data[0])
                data = re.split(", delta=", data[1])
                self.lya_alpha = float(data[0])
                data = re.split(" }", data[1])
                self.lya_delta = float(data[0]) 
            if lineStr.startswith('#linemodel solution '):
                data = re.split("for z = ", lineStr)
                data = re.split(", velocityEmission = ", data[1])
                try:
                    self.z = float(data[0]) 
                except:
                    self.z = -1
                    
                data = re.split("velocityEmission = ", lineStr)
                data = re.split(", velocityAbsorption = ", data[1])
                try:
                    self.vel_emission = float(data[0])
                    data = re.split(", merit =", data[1])
                    self.vel_absorption = float(data[0])    
                except:
                    self.vel_emission = -1
                    self.vel_absorption = -1
                
        f.close()
        self.n = len(linelambda)
        #print('len wave = {0}'.format(self.n))
        #---- default xaxis index array
        self.linelambda = range(0,self.n)
        self.linename = range(0,self.n)
        self.lineforce = range(0,self.n)
        self.linetype = range(0,self.n)  
        self.lineid = range(0,self.n)   
        self.lineamplitude = range(0,self.n) 
        self.lineerror = range(0,self.n)  
        self.lineerrorfit = range(0,self.n)  
        for x in range(0,self.n):
            self.linelambda[x] = linelambda[x]
            self.linename[x] = name[x]
            self.lineforce[x] = force[x]
            self.linetype[x] = linetype[x]
            self.lineid[x] = lineid[x]
            self.lineamplitude[x] = amplitude[x]
            self.lineerror[x] = error[x]
            self.lineerrorfit[x] = errorfit[x]
            #print("error loaded = {}".format(error[x]))
                

    def save(self, outputPath):
        text_file = open(outputPath, "w")
        data = "#\n"
        text_file.write("{}".format(data))
        
        for k in range(self.n):
            data = "{}\t{}\t{}\t-1\t-1\t{}\t-1\t-1\t\n".format(self.linetype[k], self.lineforce[k], self.linename[k], self.lineamplitude[k])
            text_file.write("{}".format(data))
        
        text_file.close()             
    
        
    def plot(self, plotType='amplitudes'):
        
        #plotValue = self.lineerror 
        if plotType=='amplitudes':
            ylabelStr = "Amplitude"
            plotValue = self.lineamplitude 
        elif plotType=='snr':
            ylabelStr = "SNR"
            plotValue = [a/self.lineerror[i] for i, a in enumerate(self.lineamplitude)]
        #
        enableYRescaleAuto = 0
        
        amax = 0.0
        emax = 0.0
        for k in range(self.n):
            if amax < plotValue[k] and self.linetype[k] == "A":
                amax = plotValue[k]
            if emax < plotValue[k] and self.linetype[k] == "E":
                emax = plotValue[k]
        aemax = max(amax, emax)
        #aemax = 2.45435e-17
        #aemax = 1e-18
        print("aemax = {}".format(aemax))
        if amax==0 and not emax==0 :
            amax = emax/20.0
        if emax==0 and not amax==0 :
            emax = amax/20.0
            
        if enableYRescaleAuto:
            plotrange = 1000.0
            coeffAmp = plotrange/aemax
            ymarginValue = 200
        else:
            coeffAmp = 1.0
            ymarginValue = aemax/5.0
            
        pp.figure(figsize=(10,8))
        linesaxisticks = []
        linesaxistickspos = []
        for k in range(self.n):
            if plotValue[k] <=0.0:
                continue
            lbl = "{} - {}A".format(self.linename[k],self.linelambda[k])
            linesaxisticks.append(lbl)
            linesaxistickspos.append(self.linelambda[k])
            x = self.linelambda[k]
            y = plotValue[k]*coeffAmp
            ccolor = 'b'
            
            ymargintext = ymarginValue+random.random()*ymarginValue*1.25
            alpha= y/emax/coeffAmp/3.0*2.0+0.3
            if self.linetype[k] == "A":
                y=-y
                ccolor = 'r'
                ymargintext = -ymargintext
                alpha= -y/amax/coeffAmp/3.0*2.0+0.33
            #print("x = {}".format(x))
            #print("y = {}".format(y))
            pp.plot((x, x), (0, y), label=self.linename[k], color=ccolor, linewidth = 3.0)
            if 1:
                pp.plot((x, x), (0, y+ymargintext), label=self.linename[k], linestyle = 'dashed', color=ccolor , alpha= alpha)
                pp.text(x, y+ymargintext, '{0}'.format(self.linename[k]), color=ccolor, alpha= alpha)
        
        
        #pp.xticks(linesaxistickspos, linesaxisticks, rotation=45,verticalalignment='top')        
        
        if enableYRescaleAuto:
            pp.ylim([-plotrange*amax/aemax - 500, plotrange*emax/aemax + 500])
        if 1:
            valf = ymarginValue*3.0
            pp.ylim([-valf, valf])
            
        pp.grid(True) # Affiche la grille
        #pp.legend(('spectrum','shifted template'), 'lower left', shadow = True)
        pp.xlabel('Rest Wavelength (Angstrom)')
        pp.ylabel('{} (x{:.3})'.format(ylabelStr, coeffAmp))
        titleStr = "{} at z={}".format(self.name, self.z)
        pp.title(titleStr) # Titre
        #pp.savefig('ExempleTrace') # sauvegarde du fichier ExempleTrace.png
        pp.show()

    def randomAmplitudes(self, coeffE=1.0, coeffA=0.5):
	#coeffE = continuum independent value: this is the 'SFR=100.0' line amplitude
        coeffE *= 1e-17 
        
        # First get OII and Halpha in range +- 30% of the SFR coefficient value
        ampOII = coeffE*(0.7+random.random()*0.6)
        ampHalpha = coeffE*(0.7+random.random()*0.6)
        
        # OIII in range +- 30% of OII
        ampOIII = ampOII*(0.7+random.random()*0.6)
        # Hbeta should be in the Balmer rule
        ampHbetaE = ampHalpha*(1.0/2.86)*(random.random())
        
        ampCIII = ampHalpha*(random.random())
        ampCIV = ampHalpha*(random.random())
	# NOTE: Emission, lineamplitude[k] = absolute amplitude of the line
        for k in range(self.n):
            if self.linetype[k] == "E":
                if self.linename[k]=="Halpha":
                    self.lineamplitude[k] = ampHalpha
                elif self.linename[k]=="Hbeta":
                    self.lineamplitude[k] = ampHbetaE
                elif self.linename[k]=="[OII]3729":
                    self.lineamplitude[k] = ampOII*0.7 
		    #coeff 0.7 = approx. doublet compensation
                elif self.linename[k]=="[OII]3726":
                    self.lineamplitude[k] = ampOII*0.7
		    #coeff 0.7 = approx. doublet compensation
                elif self.linename[k]=="[OIII](doublet-1)":
                    self.lineamplitude[k] = ampOIII
                elif self.linename[k]=="[OIII](doublet-1/3)":
                    self.lineamplitude[k] = ampOIII*0.33
                elif self.linename[k]=="LyAE":
                    #Lya in range 0-150% of OII
                    self.lineamplitude[k] = ampOII*(random.random()*1.5)
                elif self.linename[k]=="[CIII]1907":
                    #CIII in range 0-100% of Halpha
                    self.lineamplitude[k] = ampCIII
                elif self.linename[k]=="[CIII]1909":
                    self.lineamplitude[k] = ampCIII
                elif self.linename[k]=="CIV1550":
                    #CIV in range 0-100% of Halpha
                    self.lineamplitude[k] = ampCIV
                elif self.linename[k]=="CIV1548":
                    self.lineamplitude[k] = ampCIV
                elif self.linename[k]=="Hgamma":
                    #Balmer lines beyond Hgamma = 0 in emission
                    self.lineamplitude[k] = 0.0
                elif self.linename[k]=="Hdelta":
                    self.lineamplitude[k] = 0.0
                elif self.linename[k]=="Hepsilon":
                    self.lineamplitude[k] = 0.0
                elif self.linename[k]=="H8":
                    self.lineamplitude[k] = 0.0
                elif self.linename[k]=="H9":
                    self.lineamplitude[k] = 0.0
                elif self.linename[k]=="H10":
                    self.lineamplitude[k] = 0.0
                elif self.linename[k]=="H11":
                    self.lineamplitude[k] = 0.0
                else:
                    self.lineamplitude[k] = ampHalpha*(random.random()*0.33) 
                
        
        # NOTE: Absorption, lineamplitude[k] = absorption coeff of the line
        for k in range(self.n):
            if self.linetype[k] == "A":
                if self.linename[k]=="HalphaA":
                    self.lineamplitude[k] = 0.0
                elif self.linename[k]=="HbetaA":
                    #HbetaA only present if HbetaE > 0.1*HalphaE : then = in the range +-30% of the ISM coeff
                    if ampHalpha>0.0 and ampHbetaE/ampHalpha > 0.1:
                        self.lineamplitude[k] = coeffA*(0.7+random.random()*0.6)
                    else:
                        self.lineamplitude[k] = 0.0
                if self.lineforce[k] == "S":
                    #Strong ABS lines in the range 70%-130% of the ISM coeff
                    self.lineamplitude[k] = coeffA*(0.7+random.random()*0.6)
                elif self.lineforce[k] == "W":
                    #Weak ABS lines in the range 20%-30% of the ISM coeff
                    self.lineamplitude[k] = coeffA*(0.2+random.random()*0.1)
                else:
                    self.lineamplitude[k] = 0.0
              
               
        #check 
        if 0:
            for k in range(self.n):
                if self.linetype[k] == "A":               
                    self.lineamplitude[k] = 1.0
                if self.linetype[k] == "E":               
                    self.lineamplitude[k] = coeffE
         
        
    def getAmplitudeMeanNhighest(self, nhighest=5):
        ehighestpeaks = []
        ehighestpeaksIndexes = []
        for k in range(self.n):
            if self.linetype[k] == "E":
                ehighestpeaks.append(self.lineamplitude[k])
                ehighestpeaksIndexes.append(k)
                
        arr = np.array(ehighestpeaks)
        indArr = np.array(ehighestpeaksIndexes) 
        argSortedArray = arr.argsort()[::-1]
        nhighestValid = 0
        for a in indArr[argSortedArray][:nhighest]:
            if self.lineamplitude[a]>0.0:
                nhighestValid += 1
        if nhighestValid <= 0.0:
            return 0.0
        #print("highest EL indexes = {}".format(argSortedArray))
        #print("highest EL amps = {}".format(arr[argSortedArray][:nhighest]))
        highestELNames = [self.linename[a] for a in indArr[argSortedArray][:nhighestValid]]
        highestELAmps = [self.lineamplitude[a] for a in indArr[argSortedArray][:nhighestValid]]
        print("highest EL amps = {}".format(highestELAmps))
        print("highest EL names = {}".format(highestELNames))
        outVal = np.mean(highestELAmps)
        print("MEAN ( {} highest EL amps ) = {}".format(nhighestValid, outVal))
        return outVal  
        
    def getNLinesStrong(self, z=1.0, linetypefilter='E', thresSNR=3.0, thresRelAmp=3.0):
        print("\n")
        print("Linemodel-result: (N highest peaks)")
        print("thresSNR = {}".format(thresSNR))
        m_instrumentResolutionEmpiricalFactor = 230.0/325.0/2.35;
        R=250;
        validLinesIdx = []
        validLinesSnr = []
        skipnextline = 0
        for k in range(self.n): 
            if skipnextline:
                skipnextline = 0
                continue
            if self.lineerror[k]>0 and self.lineerrorfit[k]>0.0 and self.linetype[k] == linetypefilter:
                
                #check if the line doublet have higher amplitude than the individual lines        
                instrumentSigma = self.linelambda[k]*(1+z)/R*m_instrumentResolutionEmpiricalFactor;
                sigma = instrumentSigma
                print("line {}, name={}, sigma={}".format(k, self.linename[k], sigma))
                if k+1 < self.n:
                    diffLambda = abs(self.linelambda[k+1]*(1+z)-self.linelambda[k]*(1+z))
                else:
                    diffLambda = 1e32
                if k+1 < self.n and diffLambda<sigma:
                    print("name = {} : diff lambda with doublet is {} (thres={})".format(self.linename[k], diffLambda, sigma))
                    doublet_amp = self.lineamplitude[k]+self.lineamplitude[k+1]
                    doublet_error = np.mean([self.lineerror[k], self.lineerror[k+1]])
                    doublet_errorfit = np.mean([self.lineerrorfit[k], self.lineerrorfit[k+1]])
                    snr = doublet_amp/doublet_error
                    snrfit = doublet_amp/doublet_errorfit
                    if snr >= thresSNR and snrfit>=thresSNR:
                        validLinesIdx.append(k)
                        validLinesSnr.append(snr)
                        #skip the next line
                        skipnextline = 1
                else:
                    snr = self.lineamplitude[k]/self.lineerror[k]
                    snrfit = self.lineamplitude[k]/self.lineerrorfit[k]
                    print("SNR = {}, SNRFit = {}".format(snr, snrfit))
                    if snr >= thresSNR and snrfit>=thresSNR:
                        validLinesIdx.append(k)
                        validLinesSnr.append(snr)
                    
        print("Linemodel-result: (N highest peaks), found {} valid lines".format(len(validLinesIdx)))
        
        for i,k in enumerate(validLinesIdx):
            print("line {}, idx={}: lambda={}, name={}, amp={}, error={}, snr={}".format(i, k, self.linename[k], self.linelambda[k], self.lineamplitude[k], self.lineerror[k], validLinesSnr[i]))
        
        print("\n")
        
        outN = len(validLinesIdx)
        return outN
    
    def getStrongLinesCumulSNR(self, linetypeTag='E'):
        sumSnr = 0.0
        for k in range(self.n):
            if self.lineamplitude[k] > 0.0 and self.linetype[k] == linetypeTag:
                if self.lineforce[k] == 'S':
                    snr = self.lineamplitude[k]/self.lineerror[k]
                    sumSnr += snr#min(snr, 5)  
                    print("line = {}, snr = {}".format(self.linename[k], snr))
        print('Strong Lines Cumul. SNR = {}'.format(sumSnr))
        return sumSnr
        
    def getAmplitude(self, strTag, linetypeTag):
        for k in range(self.n):
            if self.linetype[k] == linetypeTag:
                if self.linename[k] == strTag:
                    outVal = self.lineamplitude[k]                
                    return outVal
        return -1
        
    def applyWeight(self, coeff):
        for k in range(self.n):
            if self.lineamplitude[k] > 0.0:
                self.lineamplitude[k] *= coeff
           
           
def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--model", dest="modelPath", default="",
                    help="path to the model to be loaded")          
          
    options = parser.parse_args()
    print(options)

    if os.path.exists(options.modelPath) :
        mpath = os.path.abspath(options.modelPath)
        print('using full path: {0}'.format(mpath))
        m = ModelResult(mpath)

        if 0:
            name_out = "201605022133_id0_z1.12688_mag20.7_sfr10.0_linemodel_amps_modified.csv"
            mpath_out = os.path.join(path,name_out)
            m.applyWeight(0.333)
            m.save(mpath_out)
    
        m.getAmplitudeMeanNhighest() 
        m.getNLinesStrong()
        m.getStrongLinesCumulSNR()
        #m.plot(plotType='amplitudes')
        m.plot(plotType='snr')
    else :
        print("Error: invalid argument count")
        exit()
    
    
    
def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()

 
if __name__ == '__main__':
    print("ModelResult")
    Main( sys.argv )
    

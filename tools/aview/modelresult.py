# -*- coding: utf-8 -*-
"""
@author: aschmitt
"""

import os
import re
import random
import matplotlib.pyplot as pp
import numpy as np

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
        for x in range(0,self.n):
            self.linelambda[x] = linelambda[x]
            self.linename[x] = name[x]
            self.lineforce[x] = force[x]
            self.linetype[x] = linetype[x]
            self.lineid[x] = lineid[x]
            self.lineamplitude[x] = amplitude[x]
            
    def plot(self):
        amax = 0.0
        emax = 0.0
        for k in range(self.n):
            if amax < self.lineamplitude[k] and self.linetype[k] == "A":
                amax = self.lineamplitude[k]
            if emax < self.lineamplitude[k] and self.linetype[k] == "E":
                emax = self.lineamplitude[k]
        aemax = max(amax, emax)
        plotrange = 1000.0
        coeffAmp = plotrange/aemax
        
        pp.figure(figsize=(10,8))
        linesaxisticks = []
        linesaxistickspos = []
        for k in range(self.n):
            if self.lineamplitude[k] <=0.0:
                continue
            lbl = "{} - {}A".format(self.linename[k],self.linelambda[k])
            linesaxisticks.append(lbl)
            linesaxistickspos.append(self.linelambda[k])
            x = self.linelambda[k]
            y = self.lineamplitude[k]*coeffAmp
            ccolor = 'b'
            ymargintext = 200+random.random()*250.0
            alpha= y/emax/coeffAmp/3.0*2.0+0.3
            if self.linetype[k] == "A":
                y=-y
                ccolor = 'r'
                ymargintext = -ymargintext
                alpha= -y/amax/coeffAmp/3.0*2.0+0.33
            #print("x = {}".format(x))
            #print("y = {}".format(y))
            pp.plot((x, x), (0, y), label=self.linename[k], color=ccolor, linewidth = 3.0)
            pp.plot((x, x), (0, y+ymargintext), label=self.linename[k], linestyle = 'dashed', color=ccolor , alpha= alpha)
            pp.text(x, y+ymargintext, '{0}'.format(self.linename[k]), color=ccolor, alpha= alpha)
        
        
        #pp.xticks(linesaxistickspos, linesaxisticks, rotation=45,verticalalignment='top')        
        
        pp.ylim([-plotrange*amax/aemax - 500, plotrange*emax/aemax + 500])
        pp.grid(True) # Affiche la grille
        #pp.legend(('spectrum','shifted template'), 'lower left', shadow = True)
        pp.xlabel('Rest Wavelength (Angstrom)')
        pp.ylabel('Amplitude (x{:.3})'.format(coeffAmp))
        pp.title(self.name) # Titre
        #pp.savefig('ExempleTrace') # sauvegarde du fichier ExempleTrace.png
        pp.show()
        
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
            
          
if __name__ == '__main__':
    path = "/home/aschmitt/code/python/simulation_perf_matrix_linemodel/amazed/output/spectrum_simu_tmp_F"
    name = "linemodelsolve.linemodel_fit_extrema_0.csv"
    mpath = os.path.join(path,name)
    print('using full path: {0}'.format(mpath))
    m = ModelResult(mpath)

    m.getAmplitudeMeanNhighest() 
    m.plot()
    
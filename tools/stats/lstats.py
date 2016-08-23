# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 19:57:47 2015

@author: aschmitt
"""
import os

import matplotlib.pyplot as pp
from matplotlib import gridspec
import numpy as np
import scipy as sp

def exportHistogram(yvect, bins, outFile="", htype="cumulative", ytype="percentage"):
    n = len(yvect)
    nbins = len(bins)
    #print 'number of bins = ' + str(nbins) + '\n'
    ybins = np.zeros(nbins)

    outsidebins = 0
    #print range(nbins-2, -1, -1)
    for x in range(0,n):
        foundbin = 0
        for ex in range(0, nbins):
            if htype=="cumulative":
                if yvect[x]<bins[ex]:
                    ybins[ex] = ybins[ex]+1
                    #break;
            elif htype =="simple":
                if ex == 0:
                    inf = bins[ex] - (bins[ex+1]-bins[ex])/2.0
                    sup = bins[ex] + (bins[ex+1]-bins[ex])/2.0
                elif ex == nbins-1:
                    inf = bins[ex] - (bins[ex]-bins[ex-1])/2.0
                    sup = bins[ex] + (bins[ex]-bins[ex-1])/2.0
                else:
                    inf = bins[ex] - (bins[ex]-bins[ex-1])/2.0
                    sup = bins[ex] + (bins[ex+1]-bins[ex])/2.0
               
                if yvect[x]>=inf and yvect[x]<sup:
                   ybins[ex] = ybins[ex]+1
                   foundbin=1
        if foundbin:
            outsidebins+=1
      
    if ytype=="percentage":
        for ex in range(0, nbins):
            ybins[ex] = ybins[ex]/float(n)*100.0
            
        
    if not outFile == "":
        # convert to percentage, and save in file 
        fout = open( outFile, "w" )  
        for ex in range(0, nbins):
            outStr = str(ex) + '\t' + str(bins[ex]) + '\t' +str(ybins[ex])
            print outStr
            fout.write( outStr  + '\n')
        strOutsideBins = "#count outside bins = {}".format(outsidebins)
        print(strOutsideBins)
        fout.write( strOutsideBins  + '\n')        
        strTotalCount = "#count total = {}".format(n)
        print(strTotalCount)
        fout.write( strTotalCount  + '\n')
        fout.close()
    
    return ybins
    
def exportHistogramComplex(yvect, mvect, bins, outFileNoExt="", htype="simple", lowhigh=[5, 25, 75, 95]):
    """
    This function calculates the histogram with the bins on mvect, 
    and then computes the average on yvect over the elements of the bin
    """
    n = len(mvect)
    if not len(yvect)==len(mvect):
        print("input data doesn't have the same length, aborting...")
    nbins = len(bins)
    #print 'number of bins = ' + str(nbins) + '\n'
    fbins = []
    for k in range(nbins):
        fbins.append([])
    #print('created empty list list: {}'.format(fbins))
    ynbins = np.zeros(nbins)
    ybinsVERYLOW = np.zeros(nbins)
    ybinsLOW = np.zeros(nbins)
    ybins49 = np.zeros(nbins)
    ybins50 = np.zeros(nbins)
    ybins51 = np.zeros(nbins)
    ybinsHIGH = np.zeros(nbins)
    ybinsVERYHIGH = np.zeros(nbins)

    #print range(nbins-2, -1, -1)
    for x in range(0,n):
        for ex in range(0, nbins):
            if htype=="cumulative":
                if mvect[x]<bins[ex]:
                    ybins[ex] = ybins[ex]+1
                    #break;
            elif htype =="simple":
                if ex == 0:
                    inf = bins[ex] - (bins[ex+1]-bins[ex])/2.0
                    sup = bins[ex] + (bins[ex+1]-bins[ex])/2.0
                elif ex == nbins-1:
                    inf = bins[ex] - (bins[ex]-bins[ex-1])/2.0
                    sup = bins[ex] + (bins[ex]-bins[ex-1])/2.0
                else:
                    inf = bins[ex] - (bins[ex]-bins[ex-1])/2.0
                    sup = bins[ex] + (bins[ex+1]-bins[ex])/2.0
               
                if mvect[x]>=inf and mvect[x]<sup:
                   fbins[ex].append(yvect[x])
                   ynbins[ex] = ynbins[ex]+1
            
    for ex in range(0, nbins):
        if len(fbins[ex])>0:
            ybinsVERYLOW[ex] = np.percentile(fbins[ex], lowhigh[0])
            ybinsLOW[ex] = np.percentile(fbins[ex], lowhigh[1])
            ybins49[ex] = np.percentile(fbins[ex], 49)
            ybins50[ex] = np.percentile(fbins[ex], 50)
            ybins51[ex] = np.percentile(fbins[ex], 51)
            ybinsHIGH[ex] = np.percentile(fbins[ex],  lowhigh[2])
            ybinsVERYHIGH[ex] = np.percentile(fbins[ex],  lowhigh[3])
        else:
            ybinsLOW[ex] = 0.0
            ybinsVERYLOW[ex] = 0.0
            ybins50[ex] = 0.0
            ybinsHIGH[ex] = 0.0
            ybinsVERYHIGH[ex] = 0.0
        
    if not outFileNoExt == "":
        # save 50 prctile in file
        fout = open( outFileNoExt+"prct50.csv", "w" )  
        for ex in range(0, nbins):
            outStr = str(ex) + '\t' + str(bins[ex]) + '\t' +str(ybins50[ex]) + '\t' +str(ynbins[ex])
            print outStr
            fout.write( outStr  + '\n')
        fout.close()  
        
        # save high prctile in file
        fout = open( outFileNoExt+"prct{}.csv".format(lowhigh[2]), "w" )  
        for ex in range(0, nbins):
            outStr = str(ex) + '\t' + str(bins[ex]) + '\t' +str(ybinsHIGH[ex]) + '\t' +str(ynbins[ex])
            print outStr
            fout.write( outStr  + '\n')
        fout.close()  
        
        # save low prctile in file
        fout = open( outFileNoExt+"prct{}.csv".format(lowhigh[1]), "w" )  
        for ex in range(0, nbins):
            outStr = str(ex) + '\t' + str(bins[ex]) + '\t' +str(ybinsLOW[ex]) + '\t' +str(ynbins[ex])
            print outStr
            fout.write( outStr  + '\n')
        fout.close()
    
    return ybins49, ybins50, ybins51, ybinsVERYLOW, ybinsLOW, ybinsHIGH, ybinsVERYHIGH, ynbins
    
    
def PlotAmazedVersusNoiseHistogram(yvect, mvect, outdir, outFileNoExt, enableExport=1, mtype='SNR'):
    if mtype=='SNR':
        print '\n\nPlotAmazedVersusNoiseHistogram:'
        vectMagBins = np.logspace(-3, 3, 50, endpoint=True)
        print 'the SNR bins are: ' + str(vectMagBins)
    else:
        print '\n\nPlotAmazedVersusMagHistogram:'
        vectMagBins = np.linspace(15, 30, 25, endpoint=True)
        print 'the mag bins are: ' + str(vectMagBins)        
    
    ybins50, ybinsVERYLOW, ybinsLOW, ybinsHIGH, ybinsVERYHIGH, ynbins = exportHistogramComplex(yvect, mvect, vectMagBins, outFileNoExt)
    enablePlot =1
    if enablePlot or enableExport:
        fig = pp.figure('Amazed stats - z error histogram versus {}'.format(mtype), figsize=(14,7))
        ax = fig.add_subplot(111)
        #pp.plot(vectMagBins, ybins25, 'b-', label="z error, prctile25")
        pp.plot(vectMagBins, ybins50, 'b-', label="z error, [5%, 25%, median, 75%, 95%]")
        #pp.plot(vectMagBins, ybins75, 'b=', label="z error, prctile75")
        pp.fill_between(vectMagBins, ybinsLOW, ybinsHIGH, alpha=0.3, label="zerror, 25%-75% prctile")
        pp.fill_between(vectMagBins, ybinsVERYLOW, ybinsVERYHIGH, alpha=0.2, label="zerror, 5%-95% prctile")
        pp.plot(vectMagBins, ynbins, 'k+', label="Number of spectra")
        #plt.semilogx()
        if mtype=='SNR':
            pp.xlim([1e-3, 1e3])
            pp.ylim([1e-5, 100])
            ax.set_xscale('log')
        else:
            pp.xlim([15, 30])
            pp.ylim([1e-5, 100])
            
        ax.set_yscale('log')
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        pp.grid(True) # Affiche la grille
        #r = Rectangle((0, 0), 1, 1) # creates rectangle patch for legend use.
        #pp.legend(r, ["z error, median", "zerror, 25%-75% prctile"],shadow=True,fancybox=True) 
        #pp.legend(shadow=True,fancybox=True, loc='center right', ncol=1)
        pp.ylabel('z error')
        pp.xlabel('{}'.format(mtype))
        name1 = "(-) z error, [Median (dark blue line), 25%-75% percentile (blue), 5%-95% percentile (light blue)]\n(+) Spectra count in bin"
        pp.title(name1)
        
        if enableExport:
            outFigFile = os.path.join(outdir,outFileNoExt + '.png')
            pp.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
        if enablePlot:
            pp.show()    
    
  
def PlotAmazedVersusBinsHistogram(yvect, mvect, outdir, outFileNoExt, enablePlot=0, enableExport=1, mtype='REDSHIFT', nPercentileDepth=2):
    if mtype=='SNR':
        print '\n\nPlotAmazedVersusNoiseHistogram:'
        vectBins = np.logspace(-3, 3, 30, endpoint=True)
        #print 'the SNR bins are: ' + str(vectBins)
    elif mtype=='MAG':
        print '\n\nPlotAmazedVersusMagHistogram:'
        #vectBins = np.linspace(15, 30, 25, endpoint=True)
        vectBins = np.linspace(15, 30, 30, endpoint=True)
        #print 'the MAG bins are: ' + str(vectBins)  
    elif mtype=='REDSHIFT':
        print '\n\nPlotAmazedVersusRedshiftHistogram:'
        #vectBins = np.linspace(0.0, 5.0, 15, endpoint=True)
        vectBins = np.linspace(0.0, 5.0, 30, endpoint=True)
        #print 'the redshift bins are: ' + str(vectBins)
    elif mtype=='SFR':
        print '\n\nPlotAmazedVersusSFRHistogram:'
        vectBins = np.logspace(-2, 2, 40, endpoint=True)
        #print 'the sfr bins are: ' + str(vectBins)  
    elif mtype=='EBMV':
        print '\n\nPlotAmazedVersusEBMVHistogram:'
        vectBins = np.logspace(-3, 2, 30, endpoint=True)
        #print 'the ebmv bins are: ' + str(vectBins)  
    elif mtype=='SIGMA':
        print '\n\nPlotAmazedVersusSIGMAHistogram:'
        vectBins = np.linspace(50, 600, 30, endpoint=True)
        #print 'the sigma bins are: ' + str(vectBins)
    elif mtype=='ERROR_SIGMA':
        print '\n\nPlotAmazedVersusERROR_SIGMAHistogram:'
        vectBins = np.linspace(-500, 500, 40, endpoint=True)
        #print 'the sigma bins are: ' + str(vectBins) 
    elif mtype=='EXT':
        print '\n\nPlotAmazedVersusEXTHistogram:'
        ext_min = np.nanmin(mvect)
        ext_max = np.nanmax(mvect)
        ext_range = ext_max-ext_min
        #ext_min -= ext_range/5.0
        #ext_min = max(ext_min, 0.0)
        #ext_max += ext_range/5.0
        if 0:
            vectBins = np.linspace(ext_min, ext_max, 40, endpoint=True)
        else:
            logmin = np.log10(ext_min)-1
            logmax = np.log10(ext_max)+1
            print("logmin = {}, logmax = {}".format(logmin, logmax))
            
            vectBins = np.logspace(logmin, logmax, 40, endpoint=True)
        #print 'the sigma bins are: ' + str(vectBins) 
    elif mtype=='LOGFHALPHA':
        print '\n\nPlotAmazedVersusLOGFHaHistogram:'
        vectBins = np.linspace(-16, -14, 20, endpoint=True)
        #print 'the sigma bins are: ' + str(vectBins)  
    
    ybins49, ybins50, ybins51, ybinsVERYLOW, ybinsLOW, ybinsHIGH, ybinsVERYHIGH, ynbins = exportHistogramComplex(yvect, mvect, vectBins, outFileNoExt)
    if enablePlot or enableExport:
        fig = pp.figure('Amazed stats - z error histogram versus {}'.format(mtype), figsize=(14,10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[12,4])
        ax2 = pp.subplot(gs[0]) #
        ax1 = pp.subplot(gs[1], sharex=ax2)  #
        
        
        nbins = len(vectBins)
        widthbins = np.zeros((nbins))
        vectBars = np.zeros((nbins))
        for i in range(nbins-1):
            widthbins[i] = (vectBins[i+1]-vectBins[i])
            vectBars[i] = (vectBins[i]-widthbins[i]/2.0)
        widthbins[nbins-1] = (vectBins[nbins-1]-vectBins[nbins-2])
        vectBars[nbins-1] = (vectBins[nbins-1]-widthbins[nbins-1]/2.0)
        #print("vectBars = {}".format(vectBars))
        #print("widthbins = {}".format(widthbins))

#        ybins50p = np.copy(ybins50);
#        ybins50m = np.copy(ybins50);
#        for i in range(nbins-1):
#            ybins50p[i] = ybins50[i]+(1.0 + 1e-6)
#            ybins50m[i] = ybins50[i]*(1.0 - 1e6)                
#        ax2.bar(vectBars, ybins50p, bottom=ybins50m, width = widthbins, label="z error, 50%")
        
        #for i in range(nbins-1):
        #    ax2.errorbar(vectBins[i], ybins50[i], 'b', xerr=0.0, yerr=widthbins[i], label="z error, 50%")
        for i in range(nbins):
            xmin = float(vectBars[i])
            xwidth = float(widthbins[i])
            xmax = xmin+xwidth
            ymin = ybins50[i]
            #ywidth = ybins50[i]*(1.0)
            ymax = ymin
            ax2.plot([xmin, xmax], [ymin, ymax], '-', color='b', fillstyle = 'full', label="z error, 50%", linewidth=2)
            #ax2.broken_barh((xmin, xwidth), (ymin, ywidth), facecolors='blue') 
        #bp = ax2.boxplot(vectBins, ybins50, patch_artist=True)
        ax2.plot(vectBins, ybins50, '+k', fillstyle = 'full', label="z error, 50%")

        data = ybinsHIGH-ybinsLOW
        mask = data.nonzero()
        #indszero = [i for i, e in enumerate(ybinsHIGH) if ybinsHIGH[i]>ybinsLOW[i]]
        #indszero = np.where(ybinsHIGH>ybinsLOW)[0]
        #print("indszero = {}".format(indszero))
        if nPercentileDepth>=1:                    
            ax2.bar(vectBars[mask], ybinsHIGH[mask], bottom=ybinsLOW[mask], width = widthbins[mask], label="z error, 25%-75%%", alpha=0.3)
            #ax2.plot(vectBins, ybinsLOW, '-b+', alpha=0.3, label="zerror, 25%-75% prctile LOW")      
            #ax2.plot(vectBins, ybinsHIGH, '-b+', alpha=0.3, label="zerror, 25%-75% prctile HIGH")
            #ax2.fill_between(vectBins, ybinsLOW, ybinsHIGH, alpha=0.3, label="zerror, 25%-75% prctile")
        
        if nPercentileDepth>=2:
            data = ybinsVERYHIGH-ybinsVERYLOW
            mask = data.nonzero()
            ax2.bar(vectBars[mask], ybinsVERYHIGH[mask], bottom=ybinsVERYLOW[mask], width = widthbins[mask], label="z error, 5%-95%", alpha=0.2)
            #ax2.plot(vectBins, ybinsVERYLOW, '-b+', alpha=0.2, label="zerror, zerror, 5%-95% prctile LOW")
            #ax2.plot(vectBins, ybinsVERYHIGH, '-b+', alpha=0.2, label="zerror, zerror, 5%-95% prctile HIGH")
            #ax2.fill_between(vectBins, ybinsVERYLOW, ybinsVERYHIGH, alpha=0.2, label="zerror, 5%-95% prctile")
        
        ax1.plot(vectBins, ynbins, 'k+', label="Number of spectra")
        d = sp.zeros(len(ynbins))        
        ax1.bar(vectBars, ynbins, bottom=np.zeros((nbins)), width = widthbins, label="count", alpha=0.4, color='k')
        #ax1.fill_between(vectBins, ynbins, where=ynbins>=d, interpolate=True, color='black', alpha=0.4)
        #ax1.fill_between(vectBins, ynbins, where=ynbins<=d, interpolate=True, color='red')

        #plt.semilogx()
        ax2.set_ylim([3e-6, 8])
        if mtype=='SNR':
            ax2.set_xlim([1e-3, 1e3])
            ax2.set_xscale('log')
            ax1.set_xscale('log')
        elif mtype=='MAG':
            #ax2.set_xlim([15, 30])
            ax2.set_xlim([15, 29])
        elif mtype=='REDSHIFT':
            ax2.set_xlim([-0.5, 5.5])
            #ax2.set_xlim([-0.5, 3.0])
        elif mtype=='SFR':
            ax2.set_xlim([1e-1, 2*1e2])
            ax2.set_xscale('log')
            ax1.set_xscale('log')
        elif mtype=='EBMV':
            ax2.set_xlim([1e-3, 1e2])
            ax2.set_xscale('log')
            ax1.set_xscale('log')
        elif mtype=='SIGMA':
            ax2.set_xlim([25, 650])
        elif mtype=='ERROR_SIGMA':
            ax2.set_xlim([-550, 550])
        elif mtype=='EXT':
            ax2.set_xlim([ext_min, ext_max])            
            ax2.set_xscale('log')
            ax1.set_xscale('log')
        elif mtype=='LOGFHALPHA':
            ax2.set_xlim([-16.5, -13.5])
        else:
            ax2.set_xlim([15, 30])
            
        ax2.set_yscale('log')
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        pp.grid(True) # Affiche la grille
        ax2.xaxis.grid(True,'major')
        ax2.yaxis.grid(True,'major')
        #r = Rectangle((0, 0), 1, 1) # creates rectangle patch for legend use.
        #pp.legend(r, ["z error, median", "zerror, 25%-75% prctile"],shadow=True,fancybox=True) 
        #pp.legend(shadow=True,fancybox=True, loc='center right', ncol=1)
        ax1.set_ylabel('COUNT', fontsize=16)
        #ax2.set_ylabel('$| z_{calc} - z_{ref} |$', fontsize=18)
        ax2.set_ylabel('| zcalc - zref |/(1+zref)', fontsize=18)
        ax1.set_xlabel('{}'.format(mtype), fontsize=18)
        name1 = "[dark blue line]=Median"
        if nPercentileDepth>=1:
             name1 = name1 + ", [blue area]=Percentile 25%-75%"
        if nPercentileDepth>=2:
             name1 = name1 + ", [light blue]=Percentile 5%-95%"
        ax2.set_title(name1)
        
        if enableExport:
            outFigFile = os.path.join(outdir,outFileNoExt + '.png')
            pp.savefig( outFigFile, bbox_inches='tight') #
        if enablePlot:
            pp.show()  
            
def PlotAmazedVersusChi2diffRatioHistogram(yvect, mvect, outdir, outFileNoExt, enableExport=1):
    print '\n\nPlotAmazedVersusBayesRatioHistogram:'
    nbins = 50
    vectBins = np.logspace(-2, 4, nbins, endpoint=True)
    #vectBins=[-1.0*a for a in reversed(vectBins)]
    #vectBins = np.linspace(-1000, 99, nbins, endpoint=True)
    print 'the bins are: ' + str(vectBins)
    
    ybins50, ybinsVERYLOW, ybinsLOW, ybinsHIGH, ybinsVERYHIGH, ynbins = exportHistogramComplex(yvect, mvect, vectBins, outFileNoExt, lowhigh=[0, 25, 75, 100])
    enablePlot =1
    if enablePlot or enableExport:
        fig = pp.figure('Amazed stats - z error histogram versus Chi2Diff.', figsize=(14,7))
        ax = fig.add_subplot(111)
        #pp.plot(vectBins, ybins25, 'b-', label="z error, prctile25")
        pp.plot(vectBins, ybins50, 'b-', label="z error, [0%, 25%, median, 75%, 100%]")
        #pp.plot(vectBins, ybins75, 'b=', label="z error, prctile75")
        pp.fill_between(vectBins, ybinsLOW, ybinsHIGH, alpha=0.3, label="zerror, 25%-75% prctile")
        pp.fill_between(vectBins, ybinsVERYLOW, ybinsVERYHIGH, alpha=0.2, label="zerror, 0%-100% prctile")
        pp.plot(vectBins, ynbins, 'k+', label="Number of spectra")
        #plt.semilogx()
        #pp.xlim([-1e3, 1e3])
        pp.ylim([1e-5, 100])
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        ax.set_xscale('log')
        ax.set_yscale('log')
        pp.grid(True) # Affiche la grille
        #r = Rectangle((0, 0), 1, 1) # creates rectangle patch for legend use.
        #pp.legend(r, ["z error, median", "zerror, 25%-75% prctile"],shadow=True,fancybox=True) 
        #pp.legend(shadow=True,fancybox=True, loc='center right', ncol=1)
        pp.ylabel('z error')
        pp.xlabel('Chi2 diff')
        name1 = "(-) z error, [Median (dark blue line), 25%-75% percentile (blue), 0%-100% percentile (light blue)]\n(+) Spectra count in bin"
        pp.title(name1)
        
        if enableExport:
            outFigFile = os.path.join(outdir,outFileNoExt + '.png')
            pp.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
        if enablePlot:
            pp.show() 
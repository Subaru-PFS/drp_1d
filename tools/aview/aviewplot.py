# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 11:34:13 2015

@author: aschmitt
"""
import os
import sys
import optparse

import matplotlib as mpl
mpl.use('Qt5Agg')
#mpl.rcParams.update(mpl.rcParamsDefault)

import matplotlib.pyplot as pp
#import seaborn as sns
#sns.set_context("paper")
##sns.set_context("poster")
#sns.set_style("whitegrid")
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


import numpy as np

import spectrum as sp
import catalog as ctlg

lines_label_method = 0 #0 = subplot with labels, #1 = using lineid_plot package on main spc view
if lines_label_method==1:
    import lineid_plot


class AViewPlot(object):
    def __init__(self, spath, npath, tpath, cpath, z, forceTplAmplitude=-1.0, forceTplDoNotRedShift = 0, enablePlot=True, scontinuumpath = ""):
        self.spath = spath
        self.scontinuumpath = scontinuumpath
        
        self.npath = npath
        self.tpath = tpath
        self.cpath = cpath
        self.z = z
        self.name = "-"
        self.spcname = "-"
        
        self.forceTplAmplitude = forceTplAmplitude;
        self.forceTplDoNotRedShift = forceTplDoNotRedShift;
        
        self.forcePlotXIndex = False
        
        # enable template plot
        if tpath=="":
            self.forcePlotNoTemplate = True
        else:
            self.forcePlotNoTemplate = False
            
        # enable noise plot
        if npath=="":
            self.forcePlotNoNoise = True
        else:
            self.forcePlotNoNoise = False
            
        # enable lines plot
        self.forcePlotNoLines = False 
        if self.cpath == "":
            self.forcePlotNoLines = True
            
        self.exportAutoDisplaysPath = ""
        
        self.load()
        if enablePlot:
            self.plot(enableReturnFig=False)
        
    def load(self):
        self.s = sp.Spectrum(self.spath) 
        
        if not self.scontinuumpath=="":
            waveErrThreshold = 0.01
            
            self.scontinuum = sp.Spectrum(self.scontinuumpath, stype='template')
            for k in range(self.scontinuum.n):
                if abs(self.s.xvect[k] - self.scontinuum.xvect[k]) > waveErrThreshold:
                    print("\nERROR: wave vector not aligned !: aborting...")
                    print("\nERROR: info: xvect[{}]={}, conti.xvect[{}]={}".format(k, self.s.xvect[k], k, self.scontinuum.xvect[k]))
                    return
                self.s.yvect[k] -=  self.scontinuum.yvect[k]
            titleStr = "spc(NoCont.)={}\n".format(self.s.name)
        else:
            titleStr = "spc={}\n".format(self.s.name)

        self.spcname = self.s.name
        
        if not self.forcePlotNoNoise:
            self.noise = sp.Spectrum(self.npath)
            titleStr += ("noise={}\n".format(self.noise.name))
        
        if not self.forcePlotNoTemplate:
            self.t = sp.Spectrum(self.tpath, stype='template')
            if self.forceTplDoNotRedShift:
                tplStr = "model"
                
                #optionally remove continuum from that model
                if True and not self.scontinuumpath=="":
                    waveErrThreshold = 0.01
                    for k in range(self.scontinuum.n):
                        if abs(self.t.xvect[k] - self.scontinuum.xvect[k]) > waveErrThreshold:
                            print("\nERROR: wave vector not aligned !: aborting...")
                            print("\nERROR: info: xvect[{}]={}, conti.xvect[{}]={}".format(k, self.t.xvect[k], k, self.scontinuum.xvect[k]))
                            return
                        self.t.yvect[k] -=  self.scontinuum.yvect[k]
                    titleStr += "(NoCont.)"
                
            else:
                tplStr = "tpl"
            titleStr += ("{}={}\n".format(tplStr, self.t.name))
                     
        titleStr += "z={}".format(self.z)
        self.name = titleStr
        
    def plot(self, enableReturnFig=False):
        #prepare plot vects  
        self.sxvect = self.s.xvect
        self.syvect = self.s.yvect
        if not self.forcePlotNoTemplate:
            if not self.forceTplDoNotRedShift:
                self.txvect = self.t.getShiftedX(self.z)
            else:
                self.txvect = self.t.xvect
            if self.forceTplAmplitude == -1.0: 
                A = self.s.ysum / self.t.ysum * 0.16
            else:
                A = self.forceTplAmplitude
            self.tyvect = self.t.getWeightedY(A)
            #prepare x,y range            
            self.xmin = self.sxvect[0]
            self.xmax = self.sxvect[self.s.n-1] 
            #self.xmin = min(self.sxvect[0], self.txvect[0])
            #self.xmax = max(self.sxvect[self.s.n-1], self.txvect[self.t.n-1])        
            self.ymin = min(min(self.syvect), min(self.tyvect))
            self.ymax = max(max(self.syvect), max(self.tyvect))
        else:
            A = 1.0
            self.xmin = self.sxvect[0]
            self.xmax = self.sxvect[self.s.n-1]       
            self.ymin = min(self.syvect)
            self.ymax = max(self.syvect)
        print("INFO: tpl Ampl. factor used is {}".format(A))
        print("INFO: tpl Ampl.*Norm [{}; {}]A is {}".format(self.xmin, self.xmax, A*self.t.getFluxNorm(self.xmin, self.xmax)))

        #compute ranges  
        if not self.exportAutoDisplaysPath == "":
            print("INFO: computing ranges for auto display")
            self.xmin = self.exportAuto_lambda_min
            self.xmax = self.exportAuto_lambda_max
            ixminSpc = self.s.getWavelengthIndex(self.xmin)
            ixmaxSpc = self.s.getWavelengthIndex(self.xmax)
            ixminTpl = self.t.getWavelengthIndex(self.xmin)
            ixmaxTpl = self.t.getWavelengthIndex(self.xmax)
            #print("exportAutoDisplaysPath: self.xmin={}".format(self.xmin))
            #print("exportAutoDisplaysPath: self.xmax={}".format(self.xmax))
            #print("exportAutoDisplaysPath: ixmin={}".format(ixmin))
            #print("exportAutoDisplaysPath: ixmax={}".format(ixmax))
            if ixmaxSpc - ixminSpc < 2:
                return
            
            if not self.forcePlotNoTemplate:
                print("computes ranges with template")
                if 0: #hybrid range
                    alpha = 0.33
                    ymintpl = min(self.tyvect[ixminTpl:ixmaxTpl])
                    yminspc = min(self.syvect[ixminSpc:ixmaxSpc])
                    self.ymin = min(ymintpl, (ymintpl+alpha*yminspc)/(1.0+alpha))
                    ymaxtpl = max(self.tyvect[ixminTpl:ixmaxTpl])
                    ymaxspc = max(self.syvect[ixminSpc:ixmaxSpc])
                    self.ymax = max(ymaxtpl, (ymaxtpl+alpha*ymaxspc)/(1.0+alpha))
                if 0: #spc+tpl range
                    self.ymin = min(min(self.syvect[ixminSpc:ixmaxSpc]), min(self.tyvect[ixminTpl:ixmaxTpl]))
                    self.ymax = max(max(self.syvect[ixminSpc:ixmaxSpc]), max(self.tyvect[ixminTpl:ixmaxTpl]))
                    
                if 0: #tpl min
                    self.ymin = min(self.tyvect[ixminTpl:ixmaxTpl])
                    self.ymax = max(self.tyvect[ixminTpl:ixmaxTpl])
                if 1: #spc min
                    self.ymin = min(self.syvect[ixminSpc:ixmaxSpc])
                    self.ymax = max(self.syvect[ixminSpc:ixmaxSpc])
            else:
                print("computes ranges without template")
                self.ymin = min(self.syvect[ixminSpc:ixmaxSpc])
                self.ymax = max(self.syvect[ixminSpc:ixmaxSpc])
                
        yrange = self.ymax-self.ymin
        xxrange = self.xmax-self.xmin
        print("Plotting aviewplot: ymin={}, ymax={}".format(self.ymin, self.ymax))

        #prepare noise vects
        if 1:
            if not self.forcePlotNoNoise:
                self.nxvect = self.noise.xvect
                self.nyvect = self.noise.yvect
                ax2Label = 'Noise'                
                useReducedAx2Range = True
                ax2Color = 'k'
        else:
            #prepare fitError vects
            tplInterp = self.t.copy()
            if not self.forceTplDoNotRedShift:
                tplInterp.applyRedshift(self.z)
            tplInterp.interpolateOnGrid(self.sxvect)
            tplInterpYvect = tplInterp.getWeightedY(A)
            self.nxvect = []
            self.nyvect = []
            for k in range(len(self.sxvect)):
                self.nxvect.append(self.sxvect[k])
                yts2 = (self.syvect[k]-tplInterpYvect[k])**2/self.noise.yvect[k]**2
                self.nyvect.append(np.sqrt(yts2))
            #self.nyvect = tplInterp.smoothGaussian(self.nyvect, 50) 
            ax2Label = 'Fit Err'       
            useReducedAx2Range = False
            ax2Color = 'b'
        
        
        
        #prepare lines        
        if not self.forcePlotNoLines:
            #self.linesx = []
            #self.linesx.append(6564)
            #self.linesx.append(5008)
            #self.linesname= []
            #self.linesname.append("Halpha")    
            #self.linesname.append("0III")   
            c = ctlg.Catalog(self.cpath)
            #print(c) 
            shiftedctlg = c.getShiftedCatalog(self.z, -1, -1) #z, "A", "S"
            #print(shiftedctlg)
            self.linesx = shiftedctlg['lambda']
            self.linestype= shiftedctlg['type']
            self.linesxrest = shiftedctlg['lambdarest']
            self.linesname = shiftedctlg['name']
            self.linesforce = shiftedctlg['force']
            #print(self.linesx)    
            #print(self.linesname)    
        
        
        #do the plotting
        if enableReturnFig:
            width = 16
            height = 15
        else:
            width = 14
            height = 11
        
        fig = pp.figure( "aview", figsize=(width, height))
        if lines_label_method==0:
            gs = gridspec.GridSpec(3, 1, height_ratios=[12,2,2])
            ax1 = pp.subplot(gs[0]) # main plotting area = spectrum
            #ax1.plot(...)
            ax2 = pp.subplot(gs[1], sharex=ax1)  # noise plotting area
            ax3 = pp.subplot(gs[2], sharex=ax1)   # lines tag area  
        elif lines_label_method == 1:
            gs = gridspec.GridSpec(2, 1, height_ratios=[12,2])
            ax1 = pp.subplot(gs[0]) # main plotting area = spectrum
            ax2 = pp.subplot(gs[1], sharex=ax1)  # noise plotting area
        
        gs.update(wspace=0.05, hspace=0.1) # set the spacing between axes. 
        # SPECTRUM
        spcColor = 'k'
        #errorColor = 'k'
        #nSigma = 12.75 #90% confidence interval bilateral
        if not self.forcePlotXIndex:
            ax1.plot(self.sxvect, self.syvect, color = spcColor, linewidth = 1, alpha=0.25, label='Flux')
            #if not self.forcePlotNoNoise:
                #error_low = np.array(self.syvect) - nSigma*np.array(self.nyvect)
                #error_high = np.array(self.syvect) + nSigma*np.array(self.nyvect)
                #ax1.fill_between(self.nxvect, error_low, error_high, color = errorColor, alpha=0.1, label="Error")
        else:
            ax1.plot(self.syvect, color = spcColor, label='Flux')
        if not self.forcePlotNoTemplate:
            if not self.forcePlotXIndex:
                ax1.plot(self.txvect, self.tyvect, 'b', label='shifted template')
            else:
                ax1.plot(self.tyvect, 'b', label='shifted template')
        #ax1.get_xaxis().set_visible(False)
        pp.setp(ax1.get_xticklabels(), visible=False)
    
        # NOISE
        if not self.forcePlotNoNoise:
            if not self.forcePlotXIndex:
                ax2.plot(self.nxvect, self.nyvect, ax2Color, label='noise')
            else:
                ax2.plot(self.nyvect, ax2Color, label='noise')

        if lines_label_method==0:
            # make ax2 nice
            pp.setp(ax2.get_xticklabels(), visible=False)
            #pp.setp(ax2.get_yticklabels(), visible=False)
            #ax2.get_yaxis().set_visible(False)
            #ax2.get_xaxis().set_visible(False)              

        #APPLY x,y range
        limmargy_up = 0.35
        limmargy_down = 0.15
        limmargx = 0.02
        ax1.set_ylim([self.ymin - yrange*limmargy_down,self.ymax+yrange*limmargy_up])
        ax1.set_xlim([self.xmin - xxrange*limmargx,self.xmax + xxrange*limmargx])
        #ax2.set_ylim([0, 10])
        #ax2.set_xlim([self.xmin,self.xmax])
        if lines_label_method==0:
            ax3.set_ylim([0, 10])
            #ax3.set_xlim([self.xmin,self.xmax])
        #ax1.ylim([self.ymin*(1-np.sign(self.ymin)*limmarg),self.ymax*(1+limmarg)])
        #ax1.xlim([self.xmin,self.xmax])             
             
        # LINES     
        line_wave = []
        line_label = []
        line_type = []
        if lines_label_method==1:
            #line_arrow_tips = []
            ak = lineid_plot.initial_annotate_kwargs()
            #print(ak)
            #ak['arrowprops']['arrowstyle'] = "->"
            pk = lineid_plot.initial_plot_kwargs()
            #print(pk)
            label_max_size = 12
        lposOffsetE = 0
        lposOffsetA = 0
        if not self.forcePlotNoLines:
            for k in range(len(self.linesx)):
                #x = self.linesx[k]*(1+self.z)
                x = self.linesx[k]
                if x<self.xmin or x>self.xmax :
                    continue
                if self.linestype[k]=='E':
                    cstyle = 'b-'
                    ccolor = 'b'
                    lpos = 5 + lposOffsetE
                    lstyle = 'dashdot'
                elif self.linestype[k]=='A':
                    cstyle = 'r-'
                    ccolor = 'r'
                    lpos = 1 + lposOffsetA
                    lstyle = 'dashed'
                if lines_label_method==0: #old method
                    ax1.plot((x, x), (-100000,100000) , cstyle, linestyle = lstyle, alpha=0.35, label=self.linesname[k] )
                    #pp.text(x, self.ymax*0.75, '{0}'.format(self.linesname[k]))
                    #print("testlinepos = {0}".format(self.t.getWeightedFlux(self.linesxrest[k], self.s.ysum / self.t.ysum)*2.0))
                    if 0 and not self.forcePlotNoTemplate:
                        ytext = self.t.getWeightedFlux(self.linesxrest[k], A)*2.0
                        #print ytext
                        try:
                            pp.text(x, ytext , '{0}'.format(self.linesname[k]))
                        except:
                            pp.text(x, self.ymax*0.75, '{0}'.format(self.linesname[k]))
                    else:
                        #showLine = self.linesforce[k]=='S'
                        showLine = True
                        if showLine:
                            #pp.text(x, self.ymax*0.75, '{0}:{1:.2f}'.format(self.linesname[k], x))
                            ax3.plot((x, x), (lpos+2, 10) , cstyle, linestyle = lstyle, alpha=0.35, label=self.linesname[k] )
                            annotation = ax3.text(x, lpos, '{0}'.format(self.linesname[k]), color=ccolor, alpha = 0.7)
                            if self.linestype[k]=='E':
                                lposOffsetE += 1
                                if lposOffsetE>3:
                                    lposOffsetE=0
                            elif self.linestype[k]=='A':
                                lposOffsetA += 1
                                if lposOffsetA>3:
                                    lposOffsetA=0
                elif lines_label_method==1: #new method with lineid_plot
                    line_wave.append(self.linesx[k])
                    
                    if len(self.linesname[k])>label_max_size:
                        labelStr = "{}...{}".format(self.linesname[k][:label_max_size], self.linesname[k][-1:])
                        #print(labelStr)
                    else:
                        labelStr = self.linesname[k]
                    line_label.append(labelStr)
                    line_type.append(self.linestype[k])
                    #if self.linestype[k]=='E':
                    #    pos = self.ymax + yrange*0.1
                    #elif self.linestype[k]=='A':
                    #    pos = self.ymin - yrange*0.25
                    #line_arrow_tips.append(pos)
                    
            if lines_label_method==1:
                
                #pk['color'] = '#f89393'
                arrow_tip_position_abs = self.ymin-yrange*self.ymax-yrange*limmargy_up/4.0
                arrow_tip_position_em = self.ymax+yrange*self.ymax+yrange*limmargy_up/4.0
                
                arrow_base_position_list_em = [f+yrange*0.075 for f in self.syvect]
                arrow_base_position_list_abs = [f-yrange*0.075 for f in self.syvect]
                print("arrow_base_position_list_em shape={}".format(len(arrow_base_position_list_em)))

                arrow_tip_position_list = []
                arrow_base_position_list = []                
                for i, l in enumerate(line_wave):
                    if 1:#line_type[i]=='E':                    
                        tip = arrow_tip_position_em
                        base = arrow_base_position_list_em[i]
                    else:
                        tip = arrow_tip_position_abs
                        base = arrow_base_position_list_abs[i]
                        
                    arrow_tip_position_list.append(tip)
                    arrow_base_position_list.append(base)
                
                print("arrow_base_position_list shape={}".format(len(arrow_base_position_list)))
                lineid_plot.plot_line_ids(self.sxvect,
                                            arrow_base_position_list_em,
                                            line_wave, line_label, 
                                            arrow_tip=arrow_tip_position_em, 
                                            annotate_kwargs=ak, plot_kwargs=pk, 
                                            ax=ax1)                
 
                for i in ax1.texts:
                    tag = i.get_label()
                    #tag = str(i.get_label())[0:-4]
                    #tag.replace('_num', '')
                    knum = tag.find("_num")
                    if knum>0:
                        #print('ktag = {}'.format(knum))
                        tag = tag[0:knum]
                    #print('tag_texts={}'.format(tag))
                    kline = line_label.index(tag)
                    if line_type[kline]=='E':
                        #ccolor = 'b'
                        ccolor = '#5f71f5' 
                    elif line_type[kline]=='A':
                        #ccolor = 'r'
                        ccolor = '#f83f3f'
                    i.set_color(ccolor)  
                    #i.set_rotation(0)
                for i in ax1.lines:
                    try:
                        tag = str(i.get_label())
                        tag = i.get_label()
                        knum = tag.find("_line")
                        if knum>0:
                            #print('ktag = {}'.format(knum))
                            tag = tag[0:knum]
                        knum = tag.find("_num")
                        if knum>0:
                            #print('ktag = {}'.format(knum))
                            tag = tag[0:knum]
                        #tag = str(i.get_label())[0:-4]
                        #tag.replace('_num', '')
                        #print('tag={}'.format(tag))
                        kline = line_label.index(tag)
                        if line_type[kline]=='E':
                            #ccolor = 'b'
                            ccolor = '#8995f5'                           
                        elif line_type[kline]=='A':
                            #ccolor = 'r'
                            ccolor = '#f89393'
                        i.set_color(ccolor)  
                        #i.set_rotation(0)
                    except:
                        continue
                        
        if lines_label_method==0:
            # make ax3 nice
            ax3.get_yaxis().set_visible(False)
               
                      
        #ax1.xaxis.set_major_locator(MultipleLocator(20))
        #ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax1.xaxis.set_minor_locator(MultipleLocator(5))
        #ax1.yaxis.set_major_locator(MultipleLocator(0.5))
        #ax1.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax1.xaxis.grid(True,'major')
        ax1.yaxis.grid(True,'major')
        
        #set the Noise axis
        if not self.forcePlotNoNoise:
            if useReducedAx2Range:
                nmean = np.mean([a for a in self.nyvect if a>0.0 and np.isfinite(a)])
                nstd = np.std([a for a in self.nyvect if a>0.0 and np.isfinite(a)])
                print("noise nmean={}, nstd={}".format(nmean, nstd))
                reducedNyvect = [a for a in self.nyvect if abs(a-nmean)<2*nstd and a>0.0]
                if len(reducedNyvect)>2:
                    nmin = min(reducedNyvect)
                    nmax = max(reducedNyvect)
                else:
                   nmin = min(self.nyvect)
                   nmax = max(self.nyvect) 
            else:
                nmin = min(self.nyvect)
                nmax = max(self.nyvect)
            #print("noise nmin={}, nmax={}".format(nmin, nmax))
            nrange = float(nmax - nmin);
            if nrange > 0.0:
                print("major locator nrange = {}".format(nrange))
                ax2.yaxis.set_major_locator(MultipleLocator(nrange/3.0))
                #ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            ax2.xaxis.grid(True,'major')
            ax2.yaxis.grid(True,'major')
            ax2.yaxis.grid(True,'minor')
            
            ax2.set_ylim([nmin, nmax])
            
        if lines_label_method==0:
            ax3.xaxis.grid(True,'major')
            ax3.yaxis.grid(False,'major')
            #pp.legend(('spectrum','shifted template'), 'lower left', shadow = True)
            if not self.forcePlotXIndex:
                ax3.set_xlabel('Angstrom')
            else:
                ax3.set_xlabel('index')
        else:
            if not self.forcePlotXIndex:
                ax2.set_xlabel('Angstrom')
            else:
                ax2.set_xlabel('index')
            
        ax1.set_ylabel('Flux')
        ax2.set_ylabel(ax2Label)
        ax1.set_title(self.name) # Titre
        
        
        if not self.exportAutoDisplaysPath == "":
            print("INFO: exporting auto display")
            if 0:
                self.xmin = self.exportAuto_lambda_min
                self.xmax = self.exportAuto_lambda_max
                ixmin = self.s.getWavelengthIndex(self.xmin)
                ixmax = self.s.getWavelengthIndex(self.xmax)
                #print("exportAutoDisplaysPath: self.xmin={}".format(self.xmin))
                #print("exportAutoDisplaysPath: self.xmax={}".format(self.xmax))
                #print("exportAutoDisplaysPath: ixmin={}".format(ixmin))
                #print("exportAutoDisplaysPath: ixmax={}".format(ixmax))
                if ixmax - ixmin < 2:
                    return
                
                if not self.forcePlotNoTemplate:
                    if 1: #hybrid range
                        alpha = 0.33
                        ymintpl = min(self.tyvect[ixmin:ixmax])
                        yminspc = min(self.syvect[ixmin:ixmax])
                        self.ymin = min(ymintpl, (ymintpl+alpha*yminspc)/(1.0+alpha))
                        ymaxtpl = max(self.tyvect[ixmin:ixmax])
                        ymaxspc = max(self.syvect[ixmin:ixmax])
                        self.ymax = max(ymaxtpl, (ymaxtpl+alpha*ymaxspc)/(1.0+alpha))
                    if 0: #spc+tpl range
                        self.ymin = min(min(self.syvect[ixmin:ixmax]), min(self.tyvect[ixmin:ixmax]))
                        self.ymax = max(max(self.syvect[ixmin:ixmax]), max(self.tyvect[ixmin:ixmax]))
                        
                    if 0: #tpl min
                        self.ymin = min(self.tyvect[ixmin:ixmax])
                        self.ymax = max(self.tyvect[ixmin:ixmax])
                    if 0: #spc min
                        self.ymin = min(self.tyvect[ixmin:ixmax])
                        self.ymax = max(self.tyvect[ixmin:ixmax])
                else:
                    self.ymin = min(self.syvect[ixmin:ixmax])
                    self.ymax = max(self.syvect[ixmin:ixmax])
                
                xxrange = self.xmax-self.xmin
                ax1.set_xlim([self.xmin - xxrange*limmargx,self.xmax + xxrange*limmargx])
                #ax2.set_xlim([self.xmin,self.xmax])
                #if lines_label_method==0:
                    #ax3.set_xlim([self.xmin,self.xmax]) 
                yrange = self.ymax-self.ymin
                ax1.set_ylim([self.ymin - yrange*limmargy_down*1.0,self.ymax+yrange*limmargy_up])
                
            outFigFile = os.path.join(self.exportAutoDisplaysPath, 'aview_{}_z{}_zoom{}.png'.format(self.spcname, self.z, self.exportAuto_name))
            #pp.savefig( outFigFile, bbox_inches='tight')
            pp.savefig( outFigFile)
            #pp.savefig('ExempleTrace') # sauvegarde du fichier ExempleTrace.png
            #pp.tight_layout()
        elif enableReturnFig:
            #gs.tight_layout(fig)
            return fig
        else:
            print("INFO: showing figure")
            pp.show()
        print '\n'
      
    def exportDisplays(self, displayExportPath):
        
        print("using z = {}".format(self.z))
        if self.z == -1 or self.z == "-1" :
            print("FOUND UNSUPPORTED REDSHIFT RESULT: Aborting export for this target...")
            return
        
        self.exportAutoDisplaysPath = displayExportPath
        zplus1 = self.z+1
        
        
        # add full range 
        smin = self.s.getWavelengthMin()/zplus1
        smax = self.s.getWavelengthMax()/zplus1
        displays_lambdas_rest_center = [(smax+smin)/2.0]
        displays_lambdas_rest_min = [smin]
        displays_lambdas_rest_max = [smax]
        displays_names = ["All"]
        
        # add lya to the pool of displayed ranges
        displays_lambdas_rest_center.append(1050)
        displays_lambdas_rest_min.append(900)
        displays_lambdas_rest_max.append(1240)
        displays_names.append("Lyb_Lyg_OVI")    
        
        # add lya to the pool of displayed ranges
        displays_lambdas_rest_center.append(1215)
        displays_lambdas_rest_min.append(1155)
        displays_lambdas_rest_max.append(1275)
        displays_names.append("Lya")   
        
        # add CIII to the pool of displayed ranges
        displays_lambdas_rest_center.append(1908)
        displays_lambdas_rest_min.append(1850)
        displays_lambdas_rest_max.append(1950)
        displays_names.append("CIII")
        
        # add CIV1550 to the pool of displayed ranges
        displays_lambdas_rest_center.append(1550)
        displays_lambdas_rest_min.append(1480)
        displays_lambdas_rest_max.append(1590)
        displays_names.append("CIV1550")
        
        # add OII_Hgamma to the pool of displayed ranges
        displays_lambdas_rest_center.append(4000)
        displays_lambdas_rest_min.append(3615)
        displays_lambdas_rest_max.append(4380)
        displays_names.append("OII_Hgamma") 
        
        # add OII to the pool of displayed ranges
        displays_lambdas_rest_center.append(4000)
        displays_lambdas_rest_min.append(3610)
        displays_lambdas_rest_max.append(3820)
        displays_names.append("OII")
        
#        # add OIII to the pool of displayed ranges
#        displays_lambdas_rest_center.append(5000)
#        displays_lambdas_rest_min.append(4925)
#        displays_lambdas_rest_max.append(5075)
#        displays_names.append("OIII")
        
        # add OIII_Hb to the pool of displayed ranges
        displays_lambdas_rest_center.append(4950)
        displays_lambdas_rest_min.append(4800)
        displays_lambdas_rest_max.append(5075)
        displays_names.append("OIII_Hb")
        
        # add Ha to the pool of displayed ranges
        displays_lambdas_rest_center.append(6562)
        displays_lambdas_rest_min.append(6450)
        displays_lambdas_rest_max.append(6780)
        displays_names.append("Ha")
        
        # add Ca to the pool of displayed ranges
        displays_lambdas_rest_center.append(3950)
        displays_lambdas_rest_min.append(3750)
        displays_lambdas_rest_max.append(4050)
        displays_names.append("CaH_CaK")
        
        ndisp = len(displays_lambdas_rest_center)
        #ndisp = 1
        for a in range(ndisp):
            self.exportAuto_lambda_center = displays_lambdas_rest_center[a]*zplus1
            self.exportAuto_lambda_min = displays_lambdas_rest_min[a]*zplus1
            self.exportAuto_lambda_max = displays_lambdas_rest_max[a]*zplus1
            self.exportAuto_name = displays_names[a]
            
            _max = max(self.s.getWavelengthMin(),self.exportAuto_lambda_min)
            _min = min(self.s.getWavelengthMax(),self.exportAuto_lambda_max)
            overlapRange = _max-_min
            print("overlapRange = {}".format(overlapRange))
            minOverlapRate = 0.66
            if overlapRange <= -(self.exportAuto_lambda_max-self.exportAuto_lambda_min)*minOverlapRate:
                print("exporting auto display for lambda rest = {}, with z = {}".format(displays_lambdas_rest_center[a], self.z))
                self.plot()
       
def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./aviewplot.py """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-s", u"--spc", help="path to the spectrum to be plotted",  dest="spcPath", default="")
    parser.add_option(u"-n", u"--noise", help="path to the noise spectrum to be plotted",  dest="noisePath", default="")
    parser.add_option(u"-t", u"--tpl", help="path to the template to be plotted",  dest="tplpath", default="")
    parser.add_option(u"-c", u"--ctlg", help="path to the catalog to be displayed",  dest="ctlgPath", default="")
    parser.add_option(u"-z", u"--redshift", help="z to be plotted",  dest="redshift", default="-1")
    parser.add_option(u"-m", u"--redshifttpl", help="force do not plot redshift",  dest="redshifttpl", default="0")
    parser.add_option(u"-a", u"--tplamp", help="amp. to be used to rescale the template",  dest="tplamp", default="-1")
    
    
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        forceTplDoNotRedShift = bool(int(options.redshifttpl))
        print("forceTplDoNotRedShift = {}".format(forceTplDoNotRedShift))
        avp = AViewPlot(options.spcPath, options.noisePath, options.tplpath, options.ctlgPath, float(options.redshift), float(options.tplamp), forceTplDoNotRedShift = forceTplDoNotRedShift )
        
        exit()
    else:
        print("Error: invalid argument count")
    
def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()  
            
if __name__ == '__main__':
    print "Catalog"
    Main( sys.argv )
        
#if __name__ == '__main__':
#    if 0:
#        path = "/home/aschmitt/data/pfs/pfs_lbg/lbgabs_1K_2z3_20J22.5"
#        name = "EZ_fits-W-F_9.fits"
#        spath = os.path.join(path,name)
#        print('using full path s: {0}'.format(spath))
#    
#        name = "EZ_fits-W-F_9.fits"
#        tpath = os.path.join(path,name)
#        print('using full path t: {0}'.format(tpath))
#        
#        avp = AViewPlot(spath, tpath, 0.1)
#    
#    if 0:
#        spath = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/1D/sc_020086957_F02P018_vmM1_red_88_1_atm_clean.fits"
#        print('Spc path is: {0}'.format(spath))
#        zval = 4.8
#        print('Redshift is: {0}'.format(zval))
#        #print(s.getRedshiftTpl(spcName))
#        tpath = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/Templates/Default/qso/SDSS_AGN.txt"
#        tpath = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/Templates/Default/galaxy/ave_Lya_no.txt"
#        print('Tpl path is: {0}'.format(tpath) )  
#        cpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/RayCatalogs/raycatalogamazedvacuum.txt"
#        print('Catalog path is: {0}'.format(cpath) )  
#        
#        #overrride zval
#        #zval = 2.11605 #for pfs lbg abs, EZ_fits-W-F_9
#        #zval = 1.1714 #for vvds1 'sc_020086397_F02P016_vmM1_red_31_1_atm_clean'
#        #zval = 1.355
#        #zval = 1.075
#        avp = AViewPlot(spath, tpath, cpath, zval)
#        
#    if 1:
#        spath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/1D/sc_020086471_F02P016_vmM1_red_107_1_atm_clean.fits"
#        spath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/1D/sc_020135682_F02P019_vmM1_red_99_1_atm_clean.fits"
#        print('Spc path is: {0}'.format(spath))
#        zval = 1.1247
#        print('Redshift is: {0}'.format(zval))
#        #print(s.getRedshiftTpl(spcName))
#        tpath = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/emission/NEW_Im_extended_blue.dat"
#        tpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/emission/NEW_Im_extended_blue.dat"
#        tpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/emission/NEW_Im_extended_blue_continuum.txt"
#        print('Tpl path is: {0}'.format(tpath) )  
#        
#        cpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/linecatalogs/linecatalogamazedair_B5.txt"
#        print('Catalog path is: {0}'.format(cpath) )  
#        
#        #overrride zval
#        #zval = 2.11605 #for pfs lbg abs, EZ_fits-W-F_9
#        #zval = 1.1714 #for vvds1 'sc_020086397_F02P016_vmM1_red_31_1_atm_clean'
#        #zval = 1.355
#        #zval = 1.075
#        avp = AViewPlot(spath, "", tpath, cpath, zval)
#    
#    if 0:
#        spath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/470026900000130.2-0.4_20_20.5_/470026900000130.2-0.4_20_20.5_EZ_fits-W-TF_0.fits"
#        print('Spc path is: {0}'.format(spath))
#        npath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/470026900000130.2-0.4_20_20.5_/470026900000130.2-0.4_20_20.5_EZ_fits-W-ErrF_0.fits"
#        print('Noise path is: {0}'.format(npath))
#        
#        zval = 0.3375
#        print('Redshift is: {0}'.format(zval))
#        #print(s.getRedshiftTpl(spcName))
#        tpath = ""
#        print('Tpl path is: {0}'.format(tpath) )  
#        cpath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/amazed/linecatalogs/linecatalogamazedvacuum_B5.txt"
#        print('Catalog path is: {0}'.format(cpath) )  
#        
#        #overrride zval
#        #zval = 2.11605 #for pfs lbg abs, EZ_fits-W-F_9
#        #zval = 1.1714 #for vvds1 'sc_020086397_F02P016_vmM1_red_31_1_atm_clean'
#        #zval = 1.355
#        #zval = 1.075
#        avp = AViewPlot(spath, npath, tpath, cpath, zval)
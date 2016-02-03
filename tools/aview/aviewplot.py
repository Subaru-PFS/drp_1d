# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 11:34:13 2015

@author: aschmitt
"""
import os

import matplotlib.pyplot as pp
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import numpy as np

import spectrum as sp
import catalog as ctlg


class AViewPlot(object):
    def __init__(self, spath, npath, tpath, cpath, z, forceTplAmplitude=-1.0, forceTplDoNotRedShift = 0, enablePlot=True):
        self.spath = spath
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
            
        self.exportAutoDisplaysPath = ""
        
        self.load()
        if enablePlot:
            self.plot()
        
    def load(self):
        self.s = sp.Spectrum(self.spath) 
        titleStr = "spc={}\n".format(self.s.name)
        self.spcname = self.s.name
        
        if not self.forcePlotNoNoise:
            self.noise = sp.Spectrum(self.npath)
            titleStr += ("noise={}\n".format(self.noise.name))
        
        if not self.forcePlotNoTemplate:
            self.t = sp.Spectrum(self.tpath, stype='template')
            if self.forceTplDoNotRedShift:
                tplStr = "model"
            else:
                tplStr = "tpl"
            titleStr += ("{}={}\n".format(tplStr, self.t.name))
                     
        titleStr += "z={}".format(self.z)
        self.name = titleStr
        
    def plot(self):
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
        print("INFO: Ampl. factor used is {}".format(A))


        #prepare noise vects
        if not self.forcePlotNoNoise:
            self.nxvect = self.noise.xvect
            self.nyvect = self.noise.yvect
        
        
        
        #prepare lines
        #self.linesx = []
        #self.linesx.append(6564)
        #self.linesx.append(5008)
        #self.linesname= []
        #self.linesname.append("Halpha")    
        #self.linesname.append("0III")   
        c = ctlg.Catalog(self.cpath)
        print(c) 
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
        fig = pp.figure( "aview", figsize=(15,11))
        gs = gridspec.GridSpec(3, 1, height_ratios=[12,2,2])
        ax1 = pp.subplot(gs[0]) # main plotting area = spectrum
        #ax1.plot(...)
        ax2 = pp.subplot(gs[1], sharex=ax1)  # noise plotting area
        ax3 = pp.subplot(gs[2], sharex=ax1)   # lines tag area     
        
        gs.update(wspace=0.05, hspace=0.1) # set the spacing between axes. 
        # SPECTRUM
        spcColor = '0.6'
        if not self.forcePlotXIndex:
            ax1.plot(self.sxvect, self.syvect, color = spcColor, label='spectrum')
        else:
            ax1.plot(self.syvect, color = spcColor, label='spectrum')
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
                ax2.plot(self.nxvect, self.nyvect, 'k', label='noise')
            else:
                ax2.plot(self.nyvect, 'k', label='noise')

        # make ax2 nice
        pp.setp(ax2.get_xticklabels(), visible=False)
        #pp.setp(ax2.get_yticklabels(), visible=False)
        #ax2.get_yaxis().set_visible(False)
        #ax2.get_xaxis().set_visible(False)
                
        # LINES     
        if not self.forcePlotNoLines:
            for k in range(len(self.linesx)):
                #x = self.linesx[k]*(1+self.z)
                x = self.linesx[k]
                if x<self.xmin or x>self.xmax :
                    continue
                if self.linestype[k]=='E':
                    cstyle = 'b-'
                    ccolor = 'b'
                    lpos = 4
                    lstyle = 'dashdot'
                elif self.linestype[k]=='A':
                    cstyle = 'r-'
                    ccolor = 'r'
                    lpos = 2
                    lstyle = 'dashed'
                ax1.plot((x, x), (-1000,1000) , cstyle, linestyle = lstyle, label=self.linesname[k] )
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
                    #pp.text(x, self.ymax*0.75, '{0}:{1:.2f}'.format(self.linesname[k], x))
                    ax3.plot((x, x), (lpos+4, 10) , cstyle, linestyle = lstyle, label=self.linesname[k] )
                    ax3.text(x, lpos, '{0}'.format(self.linesname[k]), color=ccolor)
        # make ax3 nice
        ax3.get_yaxis().set_visible(False)
               
               
        #apply x,y range
        limmarg = 0.2;
        ax1.set_ylim([self.ymin*(1-np.sign(self.ymin)*limmarg*2.0),self.ymax*(1+limmarg)])
        ax1.set_xlim([self.xmin,self.xmax])
        #ax2.set_ylim([0, 10])
        ax2.set_xlim([self.xmin,self.xmax])
        ax3.set_ylim([0, 10])
        ax3.set_xlim([self.xmin,self.xmax])
        #ax1.ylim([self.ymin*(1-np.sign(self.ymin)*limmarg),self.ymax*(1+limmarg)])
        #ax1.xlim([self.xmin,self.xmax])
        
        #ax1.xaxis.set_major_locator(MultipleLocator(20))
        #ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        #ax1.xaxis.set_minor_locator(MultipleLocator(5))
        #ax1.yaxis.set_major_locator(MultipleLocator(0.5))
        #ax1.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax1.xaxis.grid(True,'major')
        ax1.yaxis.grid(True,'major')
        
        #set the Noise axis
        if not self.forcePlotNoNoise:
            nmin = min(self.nyvect)
            nmax = max(self.nyvect)
            nrange = float(nmax - nmin);
            ax2.yaxis.set_major_locator(MultipleLocator(nrange/3.0))
            #ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            ax2.xaxis.grid(True,'major')
            ax2.yaxis.grid(True,'major')
            ax2.yaxis.grid(True,'minor')
        ax3.xaxis.grid(True,'major')
        ax3.yaxis.grid(False,'major')
        #pp.legend(('spectrum','shifted template'), 'lower left', shadow = True)
        if not self.forcePlotXIndex:
             ax3.set_xlabel('Angstrom')
        else:
            ax3.set_xlabel('index')
        ax1.set_ylabel('spectrum')
        ax2.set_ylabel('noise')
        ax1.set_title(self.name) # Titre
        
        
        if not self.exportAutoDisplaysPath == "":
            self.xmin = self.exportAuto_lambda_min
            self.xmax = self.exportAuto_lambda_max
            ixmin = self.s.getWavelengthIndex(self.xmin)
            ixmax = self.s.getWavelengthIndex(self.xmax)
            #print("ixmin={}".format(ixmin))
            #print("ixmax={}".format(ixmax))
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
            
            ax1.set_xlim([self.xmin,self.xmax])
            ax2.set_xlim([self.xmin,self.xmax])
            ax3.set_xlim([self.xmin,self.xmax])         
            ax1.set_ylim([self.ymin*(1-np.sign(self.ymin)*limmarg*10.0),self.ymax*(1+limmarg)])
            
            outFigFile = os.path.join(self.exportAutoDisplaysPath, 'aview_{}_z{}_zoom{}.png'.format(self.spcname, self.z, self.exportAuto_name))
            #pp.savefig( outFigFile, bbox_inches='tight')
            pp.savefig( outFigFile)
            #pp.savefig('ExempleTrace') # sauvegarde du fichier ExempleTrace.png
            #pp.tight_layout()
        else:
            pp.show()
        print '\n'
      
    def exportDisplays(self, displayExportPath):
        self.exportAutoDisplaysPath = displayExportPath
        zplus1 = self.z+1
        
        
        # add full range 
        smin = self.s.getWavelengthMin()/zplus1
        smax = self.s.getWavelengthMax()/zplus1
        displays_lambdas_rest_center = [(smax+smin)/2.0]
        displays_lambdas_rest_min = [smin]
        displays_lambdas_rest_max = [smax]
        displays_names = ["FullRange"]
        
        # add lya to the pool of displayed ranges
        displays_lambdas_rest_center.append(1215)
        displays_lambdas_rest_min.append(1195)
        displays_lambdas_rest_max.append(1235)
        displays_names.append("Lya")
        
        # add CIV1550 to the pool of displayed ranges
        displays_lambdas_rest_center.append(1500)
        displays_lambdas_rest_min.append(1250)
        displays_lambdas_rest_max.append(1700)
        displays_names.append("CIV1550")
        
        # add OII_Hgamma to the pool of displayed ranges
        displays_lambdas_rest_center.append(4000)
        displays_lambdas_rest_min.append(3695)
        displays_lambdas_rest_max.append(4380)
        displays_names.append("OII_Hgamma")
        
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
        displays_lambdas_rest_max.append(6650)
        displays_names.append("Ha")
        
        # add Ca to the pool of displayed ranges
        displays_lambdas_rest_center.append(3950)
        displays_lambdas_rest_min.append(3750)
        displays_lambdas_rest_max.append(4050)
        displays_names.append("CaH_CaK")
        
        
        for a in range(len(displays_lambdas_rest_center)):
            self.exportAuto_lambda_center = displays_lambdas_rest_center[a]*zplus1
            self.exportAuto_lambda_min = displays_lambdas_rest_min[a]*zplus1
            self.exportAuto_lambda_max = displays_lambdas_rest_max[a]*zplus1
            self.exportAuto_name = displays_names[a]
            
            _max = max(self.s.getWavelengthMin(),self.exportAuto_lambda_min)
            _min = min(self.s.getWavelengthMax(),self.exportAuto_lambda_max)
            overlapMin = _max-_min
            #print("overlapRange = {}".format(overlapMin))
            if overlapMin <= 2:
                print("exporting auto display for lambda rest = {}, with z = {}".format(displays_lambdas_rest_center[a], self.z))
                self.plot()
        
        
if __name__ == '__main__':
    if 0:
        path = "/home/aschmitt/data/pfs/pfs_lbg/lbgabs_1K_2z3_20J22.5"
        name = "EZ_fits-W-F_9.fits"
        spath = os.path.join(path,name)
        print('using full path s: {0}'.format(spath))
    
        name = "EZ_fits-W-F_9.fits"
        tpath = os.path.join(path,name)
        print('using full path t: {0}'.format(tpath))
        
        avp = AViewPlot(spath, tpath, 0.1)
    
    if 0:
        spath = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/1D/sc_020086957_F02P018_vmM1_red_88_1_atm_clean.fits"
        print('Spc path is: {0}'.format(spath))
        zval = 4.8
        print('Redshift is: {0}'.format(zval))
        #print(s.getRedshiftTpl(spcName))
        tpath = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/Templates/Default/qso/SDSS_AGN.txt"
        tpath = "/home/aschmitt/data/vvds/vvds2/cesam_vvds_z0_F02_DEEP/amazed/Templates/Default/galaxy/ave_Lya_no.txt"
        print('Tpl path is: {0}'.format(tpath) )  
        cpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/RayCatalogs/raycatalogamazedvacuum.txt"
        print('Catalog path is: {0}'.format(cpath) )  
        
        #overrride zval
        #zval = 2.11605 #for pfs lbg abs, EZ_fits-W-F_9
        #zval = 1.1714 #for vvds1 'sc_020086397_F02P016_vmM1_red_31_1_atm_clean'
        #zval = 1.355
        #zval = 1.075
        avp = AViewPlot(spath, tpath, cpath, zval)
        
    if 1:
        spath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/1D/sc_020086471_F02P016_vmM1_red_107_1_atm_clean.fits"
        spath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/1D/sc_020135682_F02P019_vmM1_red_99_1_atm_clean.fits"
        print('Spc path is: {0}'.format(spath))
        zval = 1.1247
        print('Redshift is: {0}'.format(zval))
        #print(s.getRedshiftTpl(spcName))
        tpath = "/home/aschmitt/data/pfs/pfs_lbg/amazed/Templates/ExtendedGalaxyEL2/emission/NEW_Im_extended_blue.dat"
        tpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/ExtendedGalaxyEL2/emission/NEW_Im_extended_blue.dat"
        tpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/Templates/linemodel/emission/NEW_Im_extended_blue_continuum.txt"
        print('Tpl path is: {0}'.format(tpath) )  
        
        cpath = "/home/aschmitt/data/vvds/vvds1/cesam_vvds_spAll_F02_1D_1426869922_SURVEY_DEEP/results_amazed/linecatalogs/linecatalogamazedair_B5.txt"
        print('Catalog path is: {0}'.format(cpath) )  
        
        #overrride zval
        #zval = 2.11605 #for pfs lbg abs, EZ_fits-W-F_9
        #zval = 1.1714 #for vvds1 'sc_020086397_F02P016_vmM1_red_31_1_atm_clean'
        #zval = 1.355
        #zval = 1.075
        avp = AViewPlot(spath, "", tpath, cpath, zval)
    
    if 0:
        spath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/470026900000130.2-0.4_20_20.5_/470026900000130.2-0.4_20_20.5_EZ_fits-W-TF_0.fits"
        print('Spc path is: {0}'.format(spath))
        npath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/470026900000130.2-0.4_20_20.5_/470026900000130.2-0.4_20_20.5_EZ_fits-W-ErrF_0.fits"
        print('Noise path is: {0}'.format(npath))
        
        zval = 0.3375
        print('Redshift is: {0}'.format(zval))
        #print(s.getRedshiftTpl(spcName))
        tpath = ""
        print('Tpl path is: {0}'.format(tpath) )  
        cpath = "/home/aschmitt/data/pfs/pfs_testsimu_20151009/amazed/linecatalogs/linecatalogamazedvacuum_B5.txt"
        print('Catalog path is: {0}'.format(cpath) )  
        
        #overrride zval
        #zval = 2.11605 #for pfs lbg abs, EZ_fits-W-F_9
        #zval = 1.1714 #for vvds1 'sc_020086397_F02P016_vmM1_red_31_1_atm_clean'
        #zval = 1.355
        #zval = 1.075
        avp = AViewPlot(spath, npath, tpath, cpath, zval)
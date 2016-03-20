# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 10:32:00 2016

@author: aschmitt
"""
from time import localtime, strftime, strptime
import matplotlib.pyplot as pp
import numpy as np
import random
import os
import PIL

class UtilBins(object):
    """
    helps exploring the parameters used for a simulation with z-mag-sfr bins
    """
    def __init__(self, target_count_per_bin):
        self.logTagStr = "UtilBins"
        self.target_count_per_bin = target_count_per_bin
        
        n_zbins_1 = 10
        n_zbins_2 = 2
        z_bins_limits_1 = np.linspace(0, 5.0, n_zbins_1+1, endpoint=True)
        z_bins_limits_2 = np.linspace(6.0, 8.0, n_zbins_2+1, endpoint=True)
        self.z_bins_limits = np.concatenate((z_bins_limits_1, z_bins_limits_2), axis=0)
        self.n_zbins = len(self.z_bins_limits)-1
        print('{}: the z bins limits are: {}'.format(self.logTagStr, str(self.z_bins_limits)))
        
        self.mag_bins_limits = [20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0]
        self.mag_bins_limits = [23.0, 24.0, 25.0, 26.0]
        self.n_mag_bins = len(self.mag_bins_limits)-1
        print('{}: the mag bins limits are: {}'.format(self.logTagStr, str(self.mag_bins_limits)))
        
        self.n_sfr_bins = 3
        self.sfr_bins_limits = np.logspace(-1.0, 2.0, self.n_sfr_bins+1, endpoint=True)
        print('{}: the sfr bins limits are: {}'.format(self.logTagStr, str(self.sfr_bins_limits)))    
        
        self.templates_path = "./templates"
        self.template_bins = []
        self.template_bins.append("NEW_Sbc_extended_extMarch2016.dat")
        self.template_bins.append("BulgedataExtensionData_extMarch2016.dat")
        self.template_bins.append("vvds_reddestdataExtensionData_extMarch2016.dat")
        self.n_template_bins = len(self.template_bins)
        print('{}: the template bins are: {}'.format(self.logTagStr, str(self.template_bins)))
    
        self.matrix = np.zeros((self.n_zbins, self.n_mag_bins, self.n_sfr_bins)) 
        
        self.datetimetag = strftime('%Y%m%d%H%M',localtime())
        self.idx_spc = 0
        
    def findIndexes(self, z, mag, sfr):
        found = False
        izFound = -1
        imagFound = -1
        isfrFound = -1
        
        for iz in range(self.n_zbins):
            if not (z > self.z_bins_limits[iz] and z <= self.z_bins_limits[iz+1]):
                continue
            
            for imag in range(self.n_mag_bins): 
                if not (mag > self.mag_bins_limits[imag] and mag <= self.mag_bins_limits[imag+1]):
                    continue
            
                for isfr in range(self.n_sfr_bins): 
                    if (sfr > self.sfr_bins_limits[isfr] and sfr <= self.sfr_bins_limits[isfr+1]):
                        izFound = iz
                        imagFound = imag
                        isfrFound = isfr
                        found = True
                        print("findIndexes, found : zbin = (i{}) [{}, {}], magbin = (i{}) [{}, {}], sfrbin = (i{}) [{}, {}]".format(izFound, self.z_bins_limits[izFound], self.z_bins_limits[izFound+1], imagFound, self.mag_bins_limits[imagFound], self.mag_bins_limits[imagFound+1], isfrFound, self.sfr_bins_limits[isfrFound], self.sfr_bins_limits[isfrFound+1]))
                        return izFound, imagFound, isfrFound
        
        return -1, -1, -1
        
    def getSpcName(self, z, mag, sfr):
        tagstr = "{}_id{}_z{:.5f}_mag{:.1f}_sfr{:.1f}".format(self.datetimetag, self.idx_spc, z, mag, sfr)
        self.idx_spc = self.idx_spc +1
        return tagstr
        
    def getPercentageFull(self):
        s = np.sum(self.matrix)
        #print s
        tot = self.target_count_per_bin * self.n_zbins * self.n_mag_bins * self.n_sfr_bins
        p = s/tot*100.0
        return p
        
    def add(self, z, mag, sfr):
        """
        this function adds an entry to the matrix if the corresponding bin is still open:
        return: 1, if the append operation was successfull, -1 if the append could not be done
        """
        iz, imag, isfr = self.findIndexes(z, mag, sfr)
        if iz==-1 and imag==-1 and isfr==-1:
            return -1        
        if self.matrix[iz, imag, isfr] < self.target_count_per_bin:
            self.matrix[iz, imag, isfr] +=1
        else:
            return -1
            
        return 1
        
    def getRandTemplate(self):
        #return self.template_bins[1]
        itpl = int(random.random()*self.n_template_bins)
        return self.template_bins[itpl]
        
    def getRandSfr(self):
        sfrmin_log = -1
        sfrmax_log = 3
        print("sfr generation between {} and {}".format(10**sfrmin_log, 10**sfrmax_log))
        sfr_log = random.random()*(sfrmax_log-sfrmin_log)+sfrmin_log
        return 10**sfr_log  
        
    def getRandVelocity(self, minV=50, maxV=500):
        print("velocity generation between {} and {}".format(minV, maxV))
        v = random.random()*(maxV-minV)+minV
        return v
        
    def getRandZMagSfr_fromOpenBin(self):
        foundOpenBin = False
        izFound = -1
        imagFound = -1
        isfrFound = -1
        
        for iz in range(self.n_zbins):
            for imag in range(self.n_mag_bins): 
                for isfr in range(self.n_sfr_bins): 
                    if self.matrix[iz, imag, isfr] < self.target_count_per_bin:
                        izFound = iz
                        imagFound = imag
                        isfrFound = isfr
                        foundOpenBin = True
                        break
                if foundOpenBin==True:
                    break
            if foundOpenBin==True:
                    break
                
        if foundOpenBin==True:
            zmin = self.z_bins_limits[izFound]
            zmax = self.z_bins_limits[izFound+1]
            z = random.random()*(zmax-zmin)+zmin
            
            magmin = self.mag_bins_limits[imagFound]
            magmax = self.mag_bins_limits[imagFound+1]
            mag = random.random()*(magmax-magmin)+magmin
            mag = 23.5
                        
            sfrmin = self.sfr_bins_limits[isfrFound]
            sfrmax = self.sfr_bins_limits[isfrFound+1]
            sfr = random.random()*(sfrmax-sfrmin)+sfrmin
            return z, mag, sfr
            
        else:
            return -1, -1, -1
            
    def exportBinCount(self, exportPath=""):
        if exportPath == "":
            print("Unbale to export : wrong path...")
        else:
            text_file = open(exportPath, "w")
            
    
    
            for iz in range(self.n_zbins):
                for imag in range(self.n_mag_bins): 
                    for isfr in range(self.n_sfr_bins):
                        count = self.matrix[iz, imag, isfr]
                        strBin = "Bin Count: zbin = (i{}) [{}, {}], magbin = (i{}) [{}, {}], sfrbin = (i{}) [{}, {}] : COUNT = {}".format(iz, self.z_bins_limits[iz], self.z_bins_limits[iz+1], imag, self.mag_bins_limits[imag], self.mag_bins_limits[imag+1], isfr, self.sfr_bins_limits[isfr], self.sfr_bins_limits[isfr+1], count)
                        #print("{}".format(strBin)) 
                        text_file.write("{}\n".format(strBin))
                    
            text_file.close()
        
    def exportBinCountPlot(self, outputDirectory=""):
        if outputDirectory == "":
            print("Unbale to export : wrong path...")
            return
            
        for iz in range(self.n_zbins):
            countmatrix = self.matrix[iz]
            cmap = pp.get_cmap('Greys')
            cmin = 0
            cmax = self.target_count_per_bin
        
            pp.clf()
            pp.close()
            pp.figure('bin count')
            pp.imshow(np.transpose(countmatrix), interpolation='nearest',aspect='auto', cmap=cmap, alpha=0.8)
            pp.clim(cmin,cmax)
            
            pp.xlabel('MAG')
            pp.ylabel('SFR')
            mag_axis_ticks = ["{:.1f}".format(a) for a in self.mag_bins_limits[::]]
            mag_axis_ticks_inds = [int(a)-0.5 for a in range(self.n_mag_bins)[::]]
            pp.xticks(mag_axis_ticks_inds, mag_axis_ticks)
            
            sfr_axis_ticks = ["{:.1f}".format(a) for a in self.sfr_bins_limits[::]]
            sfr_axis_ticks_inds = [int(a)-0.5 for a in range(self.n_sfr_bins)[::]]
            pp.yticks(sfr_axis_ticks_inds, sfr_axis_ticks)
        
            zmin = self.z_bins_limits[iz]
            zmax = self.z_bins_limits[iz+1]
            name1 = "Count matrix for {} spectra in z range [ {:.1f} ; {:.1f} ]".format(self.target_count_per_bin, zmin, zmax)
            pp.title(name1)
            pp.grid(True,'major',linewidth=1)
            pp.colorbar()
        
            outFileNoExt = "count_zmin{:.1f}_zmax{:.1f}".format(zmin, zmax)
            outFigFile = os.path.join(outputDirectory,outFileNoExt + '.png')
            pp.savefig( outFigFile, bbox_inches='tight')
            #pp.show()
        
        mergeImagesInFolder(outputDirectory, ncols = 5)
        
def mergeImagesInFolder(dir_path, ncols):
    merge_image_name = "bincount.png"    
    log_tag = "Image merge: "
    print("\n\nINFO: Now merging performance images...")
    image_files = []
    for file in sorted(os.listdir(dir_path)):
        if file.endswith(".png") and not str(file)==merge_image_name:
            image_files.append(os.path.join(dir_path, file))
            print("{}found image {}".format(log_tag, file))
        
    n_images = len(image_files)
    print("{}file list N = {}".format(log_tag, n_images))

    if len(image_files)<1:
        print("{}No Images Found ... aborting".format(log_tag))
        return
        
    #prepare big merged image         
    nx = ncols
    if n_images<ncols:
        nx = n_images
    ny = int(n_images//ncols)+int(n_images%ncols>0)
    #print("q = {}\nr = {}".format( int(n_images//ncols), int(n_images%ncols>0)))

    image = PIL.Image.open(image_files[0])
    lx = image.size[0]
    ly = image.size[1]
    merged_image = PIL.Image.new("RGB", (nx*lx, ny*ly), "white")
    
    ix = 0
    iy = 0
    for f in image_files:
        image = PIL.Image.open(f)
        #print image.size
        merged_image.paste(im=image, box=(ix, iy))
        ix = ix+lx
        if ix>=ncols*lx:
            iy = iy+ly
            ix = 0
    
    merged_image.save(os.path.join(dir_path, merge_image_name))
            
            

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 

@author: aschmitt
"""
import os
import matplotlib as mpl
import matplotlib.pyplot as pp
from matplotlib import gridspec
import matplotlib.cm as cm
import numpy as np
import scipy as sp
  
import PIL
  
def exportPerformances(dataDiff, outputDirectory):

    # dataDiff expected cols are:
    #Spectrum ID	MAGI	ZREF	ZFLAG	ZCALC	MERIT	TPL	METHOD    SNR	SFR E(B-V) Sigma DIFF
    izref = 2

    print '\n\nexportPerformances:'
    n_zbins_1 = 10
    n_zbins_2 = 2
    z_bins_limits_1 = np.linspace(0, 5.0, n_zbins_1+1, endpoint=True)
    z_bins_limits_2 = np.linspace(6.0, 8.0, n_zbins_2+1, endpoint=True)
    z_bins_limits = np.concatenate((z_bins_limits_1, z_bins_limits_2), axis=0)
    n_zbins = len(z_bins_limits)-1
    print 'the z bins limits are: ' + str(z_bins_limits)  
    

    n = len(dataDiff)
    print("datadiff n = {}".format(n))
    for izbin in range(len(z_bins_limits)-1):
        zbin_dataset = []
        for x in range(n):
            zref = dataDiff[x][izref]
            if zref >= z_bins_limits[izbin] and zref < z_bins_limits[izbin+1]:
                zbin_dataset.append(dataDiff[x])
                #print("line kept is {}".format(dataDiff[x]))
        print("\n\nThis z bin n = {}".format(len(zbin_dataset)))
        plotMagSFRPerformanceMatrix(zbin_dataset, z_bins_limits[izbin], z_bins_limits[izbin+1], 1e-3, 0.8, outputDirectory)
    mergeImagesInFolder(outputDirectory, ncols = 4)
      
def mergeImagesInFolder(dir_path, ncols):
    merge_image_name = "performance.png"    
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
     
     
def plotMagSFRPerformanceMatrix(bin_dataset, zmin, zmax, catastrophic_failure_threshold, success_rate_thres, outdir, enableExport=True):    
    """
    bin_dataset : expected cols are  [ #Spectrum ID	MAGI	ZREF	ZFLAG	ZCALC	MERIT	TPL	METHOD    SNR	SFR E(B-V) Sigma DIFF ]
    success_rate_thres : should be in range [0; 1]
    """
    #n_mag_bins = 14
    #mag_bins_limits = np.linspace(16.0, 30.0, n_mag_bins+1, endpoint=True)
    mag_bins_limits = [12.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 32.0]
    n_mag_bins = len(mag_bins_limits)-1
    print 'the mag bins limits are: ' + str(mag_bins_limits)  
    imag = 1
    mag_axis_ticks = ["{:.1f}".format(a) for a in mag_bins_limits[::]]
    mag_axis_ticks_inds = [int(a)-0.5 for a in range(n_mag_bins)[::]]
    
    n_sfr_bins = 4
    sfr_bins_limits = np.logspace(-1.0, 3.0, n_sfr_bins+1, endpoint=True)
    print 'the sfr bins limits are: ' + str(sfr_bins_limits) 
    isfr = 9
    sfr_axis_ticks = ["{:.1f}".format(a) for a in sfr_bins_limits[::]]
    sfr_axis_ticks_inds = [int(a)-0.5 for a in range(n_sfr_bins)[::]]
    
    n = len(bin_dataset)                
    matrix = np.zeros((n_mag_bins,n_sfr_bins))
    n_matrix = np.zeros((n_mag_bins,n_sfr_bins))
    success_rate_matrix = np.zeros((n_mag_bins,n_sfr_bins)) 
    for imbin in range(n_mag_bins):
        m_bin_dataset = []
        for x in range(n):
            mag = bin_dataset[x][imag]
            #print("imbin = {}".format(imbin))
            #print("mag = {}".format(mag))
            #print("mag bin min = {}".format(mag_bins_limits[imbin]))
            #print("mag bin max = {}".format(mag_bins_limits[imbin+1]))
            if mag >= mag_bins_limits[imbin] and mag < mag_bins_limits[imbin+1]:
                m_bin_dataset.append(bin_dataset[x])
                #print("line kept is {}".format(bin_dataset[x]))
                
        n_mag = len(m_bin_dataset)
        for isfrbin in range(n_sfr_bins):
            sfr_bin_dataset = []
            for x in range(n_mag):
                sfr = m_bin_dataset[x][isfr]
                if sfr >= sfr_bins_limits[isfrbin] and sfr < sfr_bins_limits[isfrbin+1]:
                    sfr_bin_dataset.append(m_bin_dataset[x])
                 
            n_matrix[imbin][isfrbin] = len(sfr_bin_dataset) 
            success_rate_matrix[imbin][isfrbin] = getSuccessRate(sfr_bin_dataset, catastrophic_failure_threshold)            
            
            if len(sfr_bin_dataset) == 0:
                matrix[imbin][isfrbin] = 0
            elif success_rate_matrix[imbin][isfrbin] > success_rate_thres:
                matrix[imbin][isfrbin] = -1
            else:
                matrix[imbin][isfrbin] = 1
            
    pp.clf()
    pp.close()
    pp.figure('performance')
    #fig = pp.figure('catalog z map', figsize=(9, 8))
    #ax = fig.add_subplot(111)
    #i = ax.matshow(np.transpose(matrix), interpolation='nearest', aspect='equal', cmap=cmap)
    
    cmap = pp.get_cmap('seismic')
    cmin = -1.4
    cmax = 1.4
    
    use_transparency = True
    if use_transparency: 
        max_alpha_map = 75.0
        alpha_map = np.zeros((n_mag_bins,n_sfr_bins))
        for imbin in range(n_mag_bins):
            for isfrbin in range(n_sfr_bins):
                _count = n_matrix[imbin][isfrbin]
                if not _count == 0:
                    _modified_count = np.log(_count+1)*20
                    _modified_count = min(max_alpha_map, float(_modified_count))
                else:
                    _modified_count = 0
                alpha_map[imbin][isfrbin] = _modified_count/float(max_alpha_map)/1.25 #generate transparency
                
        #generates rgbA matrix for transparency
        norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        tmp = m.to_rgba(np.transpose(matrix))
        print("tmp shape is {}".format(np.shape(tmp)))
        
        for imbin in range(n_mag_bins):
            for isfrbin in range(n_sfr_bins):
                alpha = alpha_map[imbin,isfrbin]
                #tmp[imbin,isfrbin] = cmap(matrix[imbin,isfrbin])
                #tmp[isfrbin,imbin][3] = alpha
                for c in range(3):
                    tmp[isfrbin,imbin][c] = (alpha)*tmp[isfrbin,imbin][c]+(1-alpha)
               
        #pp.imshow(np.transpose(matrix), interpolation='nearest',aspect='auto', cmap=cmap, alpha=0.8)
        pp.imshow(tmp, aspect='auto', interpolation='none', alpha=1.0)
    else:
        pp.imshow(np.transpose(matrix), interpolation='nearest',aspect='auto', cmap=cmap, alpha=0.8)
        pp.clim(cmin,cmax)
    
    pp.xticks(mag_axis_ticks_inds, mag_axis_ticks)#, rotation='horizontal',verticalalignment='bottom')
    pp.yticks(sfr_axis_ticks_inds, sfr_axis_ticks)#, rotation='horizontal',verticalalignment='bottom')
    
    for imbin in range(n_mag_bins):
        for isfrbin in range(n_sfr_bins):
            count = int(n_matrix[imbin][isfrbin])
            if count >0:
                c = "{}".format(count) 
                color_change_threshold = 0.7
                if alpha_map[imbin,isfrbin] < color_change_threshold:
                    color_per = 'darkred'
                else:
                    color_per = 'lightpink'
                
                if success_rate_matrix[imbin][isfrbin] > success_rate_thres:
                    if alpha_map[imbin,isfrbin] < color_change_threshold:
                        color_per = 'darkblue'
                    else:
                        color_per = 'lightskyblue'
                #print("drawing text {} at ({},{})".format(c, mag_axis_ticks_inds[imbin], sfr_axis_ticks_inds[isfrbin]))
                if alpha_map[imbin,isfrbin] < color_change_threshold:
                    color_count = 'black'
                else:
                    color_count = 'white'
                pp.text(mag_axis_ticks_inds[imbin]+0.5, sfr_axis_ticks_inds[isfrbin]+0.5, c, va='center', ha='center',fontsize=14, color=color_count)
                per = "{}%".format(int(100*success_rate_matrix[imbin][isfrbin]+0.5)) 
                #print("drawing text {} at ({},{})".format(c, mag_axis_ticks_inds[imbin], sfr_axis_ticks_inds[isfrbin]))
                pp.text(mag_axis_ticks_inds[imbin]+0.5, sfr_axis_ticks_inds[isfrbin]+0.7, per, va='center', ha='center',fontsize=12, color=color_per)
    
    #pp.subplots_adjust(bottom=0.15)
    
    pp.xlabel('MAG')
    pp.ylabel('SFR')
    name1 = "Performance matrix for {} spectra in z range [ {:.1f} ; {:.1f} ]\nFailure threshold dz/z > {}\nSuccess(Blue) strict(>) threshold: {}%".format(n, zmin, zmax, catastrophic_failure_threshold, int(success_rate_thres*100))
    pp.title(name1)
    
    #pp.colorbar()
    #ml = MultipleLocator(5)
    #pp.axes().yaxis.set_minor_locator(ml)

    pp.grid(True,'major',linewidth=1)
    #pp.grid(True,'minor', linewidth=1)
    if enableExport:
        outFileNoExt = "performance_zmin{:.1f}_zmax{:.1f}".format(zmin, zmax)
        outFigFile = os.path.join(outdir,outFileNoExt + '.png')
        pp.savefig( outFigFile, bbox_inches='tight')
    else:
        pp.show() 
  
      
def getSuccessRate(datasetdiff, catastrophic_failure_threshold):
    idiff = 12
    n_total = len(datasetdiff)
    if n_total == 0:
        return 0.0
    n_success = 0
    for e in datasetdiff:
        if abs(e[idiff]) < catastrophic_failure_threshold:
            n_success += 1
    
    rate = float(n_success)/float(n_total)
    return rate
    
# v3.0    2016-07-06 added stats support for euclid-simu dataset type.
# v2.1	2015-07-06 stats hist saved (txt export) in 2 different bin sizes
# v2.02	2015-06-19 stats histograms different colors with methods/templates subsets
# v2.01	2015-06-18 plot different colors for methods/template subsets
# v2.0	2015-06-18 allow calc file to be in a different order than the ref file. reordering is done using the spectrum name 
# v1.52 2015-04-15 add output file failures.txt listing the failures over a threshold (default 0.01)
# v1.51 2015-04-13 improved robustness for PFS data loading and new amazed processing
# v1.5	2015-04-10 creation of processAmazedOutputVVDS.py (replacement of original processAmazedOutputVVDS.py)
# v1.21 2015-04-10 set stats files output path to calculated input file path
# v1.2 	2015-04-09 diff graphs normalized by 1+z
# v1.1 	2015-03-31
#!/usr/bin/python

import sys
import os
import optparse

from astropy.io import ascii
import numpy as np

from astropy.io import fits
import random
import matplotlib.pyplot as pp
import matplotlib.cm as cm
from matplotlib.patches import Rectangle

import lstats
import lperf

# global variables, default for VVDS
iRefID = 0
iRefZ = 4
iRefMag = 6
iRefFlag = 5
iRefSFR = -1
iRefEBmV = -1
iRefSigma = -1
iRefLogHalpha = -1

def setMuseRefFileType():
    global iRefZ, iRefMag, iRefFlag, iRefSFR, iRefEBmV, iRefSigma
    iRefZ = 1
    iRefMag = -1
    iRefFlag = -1
    
    iRefSFR = -1
    iRefEBmV = -1
    iRefSigma = -1
    
def setVUDSRefFileType():
    global iRefZ, iRefMag, iRefFlag
    iRefZ = 3
    iRefMag = 5
    iRefFlag = 4  
    
def setVVDSRefFileType():
    global iRefZ, iRefMag, iRefFlag
    iRefZ = 4
    iRefMag = 6
    iRefFlag = 5

def setVVDS2RefFileType():
    global iRefZ, iRefMag, iRefFlag
    iRefZ = 5
    iRefMag = 4
    iRefFlag = 6

def setPFSRefFileType():
    global iRefZ, iRefMag, iRefFlag, iRefSFR, iRefEBmV, iRefSigma
    iRefZ = 1
    iRefMag = 2
    iRefFlag = -1
    iRefSFR = 5
    iRefEBmV = 4
    iRefSigma = 6
    
def setKeckRefFileType():
    global iRefID, iRefZ, iRefMag, iRefFlag, iRefSFR, iRefEBmV, iRefSigma
    iRefID = 3
    iRefZ = 1
    iRefMag = 2
    iRefFlag = -1
    iRefSFR = 5
    iRefEBmV = 4
    iRefSigma = 6

def setSIMULMRefFileType():
    global iRefZ, iRefMag, iRefFlag, iRefSFR, iRefEBmV, iRefSigma
    iRefZ = 1
    iRefMag = 2
    iRefFlag = -1
    iRefSFR = 5
    iRefEBmV = 4
    iRefISM = 6
    iRefSigma = 7
    
def setSIMUEuclid2016RefFileType():
    global iRefZ, iRefMag, iRefFlag, iRefSFR, iRefEBmV, iRefSigma, iRefLogHalpha
    iRefZ = 1
    iRefMag = 2
    iRefFlag = -1
    iRefSFR = 5
    iRefEBmV = 4
    iRefISM = -1
    iRefSigma = 6
    iRefLogHalpha = 7
    

def ProcessDiff( refFile, calcFile, outFile, reftype ) :
    global iRefID, iRefZ, iRefMag, iRefFlag, iRefSFR, iRefEBmV, iRefSigma, iRefLogHalpha   
     

    #*************** open ref file
    if 0:     
        fref = open(refFile, 'r')
        dataRefStr = fref.read()
        fref.close()
        dataRef = ascii.read(dataRefStr)
    else:
        dataRef = loadRef( refFile, reftype );
    dataRef_names = [a[iRefID] for a in dataRef]
    print("ProcessDiff : Dataref N = {}".format(len(dataRef)))    
    print("ProcessDiff : Dataref, first elt: {}".format(dataRef[0]))
    #print("ProcessDiff : Dataref, second elt: {}\n".format(dataRef[1]))

    #*************** open calc file
    if 0: 
        fcalc = open(calcFile, 'r')
        dataCalcStr = fcalc.read()
        fcalc.close()
        #dataCalc = ascii.read(dataCalcStr)
        #dataCalc = ascii.read(calcFile, data_start=1, delimiter='\t');
        #dataCalc = np.loadtxt(calcFile, delimiter='\t', usecols=(1, 3), unpack=True)
        dataCalc_raw = ascii.read(dataCalcStr, data_start=0, delimiter="\t")
        #print dataCalc[0]
    else:
        dataCalc_raw = loadCalc( calcFile );
    dataCalc_names_raw = [a[0] for a in dataCalc_raw]
    
    ## PFS specific: filter dataCalc_names_raw to exclude ambiguity due to input fits names
    for i,x in enumerate(dataCalc_names_raw):
        excludeStr = "FILTERED_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version"
        if x.find(excludeStr)==-1:
            excludeStr = "FILTERED_FEB2015_1stapprox_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version"
        if x.find(excludeStr)==-1:
            excludeStr = "FILTERED_SEPT2015_1loopw100_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version"
        if x.find(excludeStr)==-1:
            excludeStr = "FILTERED_SEPT2015_2loopw50_SET_reallyjustlinecont_1k_0.5z1.8_0.1A100_Version"
        if x.find(excludeStr)==-1:
            excludeStr = "470026900000130.2-0.4_20_20.5_EZ_fits-W-TF"
        if x.find(excludeStr)==-1:
            excludeStr = "470026900000130.2-0.4_20_20.5_EZ_fits-W-F"
        if x.find(excludeStr)==-1:
            excludeStr = "SPC_fits-W-F_"
        if x.find(excludeStr)==-1:
            excludeStr = "SPC_fits-W-TF_"
        if x.find(excludeStr)==-1:
            excludeStr = "SPC_fits-W-FILT_"
        
        dataCalc_names_raw[i] = x.replace(excludeStr, "")
    ##
        
    
    print("ProcessDiff : dataCalc_raw N = {}".format(len(dataCalc_raw)))    
    print("ProcessDiff : dataCalc_raw, first elt: {}".format(dataCalc_raw[0]))
    #print("ProcessDiff : dataCalc_raw, second elt: {}\n".format(dataCalc_raw[1]))
      
    #*************** reorder calc data by filename
    inds = []
    #readonlyN = 300;
    for s in dataRef_names:#[0:readonlyN]:
        #remove extension .fits from pandora results
        #s = s[0:-5]
        #p = [i for i,x in enumerate(dataCalc_names_raw) if str(s) == x] # PFS batch 1 to 5
        p = [i for i,x in enumerate(dataCalc_names_raw) if str(s) in x] # PFS batch 6 onwards
        if len(p) == 1:
            #print "OK : index found : ref={0} and calc={1}".format(s,p)
            inds.append(p[0])
        elif len(p) > 1:
            print "ERROR : multiple calc index found for a given ref name!!".format()
            print "ERROR : dataRef_names = {}\n".format(s)
            print "ERROR : dataCalc_names_raw = {}\n".format(p)
            stop
        else:
            inds.append(-1)
            print "ERROR : index not found while comparing ref and calc lists!!".format()
            print "ERROR : dataRef_names = {}\n".format(s)
            stop
    #print inds
    dataCalc = [dataCalc_raw[i] for i in inds]
    # end re-order
    
    ### try to open snr file
    snrFile = os.path.join(os.path.dirname(os.path.abspath(refFile)), "snr2_TF_ErrF.csv")
    print("SNRfile = {}".format(snrFile))
    if 1 and os.path.exists(snrFile):
        fsnr = open(snrFile, 'r')
        dataSnrStr = fsnr.read()
        fsnr.close()
        dataSnr_raw = ascii.read(dataSnrStr)
        dataSnr_names = [a[0].replace("SPC_fits-W-TF_", "").replace(".fits", "") for a in dataSnr_raw]
        print("dataSnr, first elt: {}".format(dataSnr_raw[0]))
        #*************** reorder snr data by filename
        inds = []
        #readonlyN = 300;
        for s in dataRef_names:#[0:readonlyN]:
            #remove extension .fits from pandora results
            #s = s[0:-5]
            print("dataRef_names entry: {}".format(s))
            print("dataSnr_names first entry: {}".format(dataSnr_names[0]))
        
            #p = [i for i,x in enumerate(dataSnr_names) if str(s) == x] # PFS batch 1 to 5
            p = [i for i,x in enumerate(dataSnr_names) if str(s) in x]# PFS batch 6 onwards
            if len(p) >0:
                #print("p={}".format(p))
                #print("x={}".format(dataSnr_names[p[0]]))
                print "OK : index found : ref={0} and snr={1}".format(s,dataSnr_names[p[0]])
                inds.append(p[0])
            else:
                inds.append(-1)
                print "ERROR : index not found : {0}".format(s)
                stop
        #print inds
        dataSnr = [dataSnr_raw[i] for i in inds]
        # end re-order
    else:
        dataSnr = []
        
    ### try to open external data file :it has to be placed in the ref file directory
    externalFile = os.path.join(os.path.dirname(os.path.abspath(refFile)), "external.csv")
    print("Externalfile = {}".format(externalFile))
    if 1 and os.path.exists(externalFile):
        fext = open(externalFile, 'r')
        dataExtStr = fext.read()
        fext.close()
        dataExt_raw = ascii.read(dataExtStr)
        dataExt_names = [a[0].replace("SPC_fits-W-TF_", "").replace(".fits", "") for a in dataExt_raw]
        print("dataExt, first elt: {}".format(dataExt_raw[0]))
        #*************** reorder ext data by filename
        inds = []
        #readonlyN = 300;
        for s in dataRef_names:#[0:readonlyN]:
            #remove extension .fits from pandora results
            #s = s[0:-5]
            print("dataRef_names entry: {}".format(s))
            print("dataExt_names first entry: {}".format(dataExt_names[0]))
        
            #p = [i for i,x in enumerate(dataSnr_names) if str(s) == x] # PFS batch 1 to 5
            p = [i for i,x in enumerate(dataExt_names) if str(s) in x]# PFS batch 6 onwards
            if len(p) == 1:
                #print("p={}".format(p))
                #print("x={}".format(dataSnr_names[p[0]]))
                print "OK : index found : ref={0} and snr={1}".format(s,dataExt_names[p[0]])
                inds.append(p[0])
            else:
                inds.append(-1)
                print "ERROR : index not found : {0}".format(s)
                stop
        #print inds
        dataExt = [dataExt_raw[i] for i in inds]
        # end re-order
    else:
        print('EXT file not found')
        dataExt = []
    
    

    print "INFO: ref file size = " + str(len(dataRef))
    print "INFO: calc file size = " + str(len(dataCalc))
    n = min(len(dataRef), len(dataCalc))
    print "INFO: n set to = " + str(n)
            
    #Write the diff file
    f = open( outFile, "w" )
    f.write( "#ID\tMAGI\tZREF\tZFLAG\tZCALC\tMERIT\tTPL\tMETHOD\tSNR\tSFR\tE(B-V)\tSigma\tLogHalpha\tExtValue\tDIFF\n" )
    for k in range(0,n):
        if iRefFlag>-1:
            flagValStr = str(dataRef[k][iRefFlag])
        else:
            flagValStr = "-1"
        if iRefSFR>-1:
            SFRValStr = str(dataRef[k][iRefSFR])
        else:
            SFRValStr = "-1"
        if iRefEBmV>-1:
            EBmVValStr = str(dataRef[k][iRefEBmV])
        else:
            EBmVValStr = "-1"
        if iRefSigma>-1:
            SigmaValStr = str(dataRef[k][iRefSigma])
        else:
            SigmaValStr = "-1"
        if iRefLogHalpha>-1:
            LogHalphaValStr = str(dataRef[k][iRefLogHalpha])
        else:
            LogHalphaValStr = "-1"
            
        if len(dataSnr)>0:
            snrValStr = str(dataSnr[k][1])
        else:
            snrValStr = "-1" 
            
        if len(dataExt)>0:
            extValStr = str(dataExt[k][1])
        else:
            extValStr = "-1"       

        tplStr = "-1" #str(dataCalc[k][4])
        zdiff = dataCalc[k][1]-dataRef[k][iRefZ]
                
        
        if 0:
            print "\n"
            print("DatarefNames: {}".format(dataRef_names[k]))
            print("Name: {}".format(str(dataCalc[k][0])))
            print("zref: {}".format(str(dataRef[k][iRefZ]) ))
            print str(dataRef[k][iRefMag])
            
            print("zcalc: {}".format(str(dataCalc[k][1])))
            print str(dataCalc[k][2])
            print str(dataCalc[k][3])
	        #print str(dataCalc[k][4])
            print("diff: {}".format(str(dataRef[k][iRefZ] - dataCalc[k][1])))
            print("snr: {}".format(snrValStr))
            print("ext: {}".format(extValStr))
            
         
        #f.write( str(dataCalc[k][0]) + "\t" + str(dataRef[k][iRefMag]) + "\t" + str(dataRef[k][iRefZ]) + "\t" + str(dataRef[k][iRefFlag]) + "\t" + str(dataCalc[k][1]) + "\t" + str(dataCalc[k][2]) + "\t" + str(dataCalc[k][3]) + "\t" + str(dataCalc[k][4]) + "\t" + str(dataRef[k][iRefZ] - dataCalc[k][1]) + "\n" )
        f.write( str(dataCalc[k][0]) + "\t")
        f.write( str(dataRef[k][iRefMag]) + "\t")
        f.write( str(dataRef[k][iRefZ]) + "\t")
        f.write( flagValStr + "\t")
        f.write( str(dataCalc[k][1]) + "\t")
        f.write( str(dataCalc[k][2]) + "\t")
        f.write( str(dataCalc[k][3]) + "\t")
        f.write( tplStr + "\t")
        f.write( snrValStr + "\t")
        f.write( SFRValStr + "\t")
        f.write( EBmVValStr + "\t")
        f.write( SigmaValStr + "\t")
        f.write( LogHalphaValStr + "\t")
        f.write( extValStr + "\t" )
        f.write( str(zdiff) + "\n" )
    
    f.close()

    if 0:
        #*************** create subsets
        subset_index = 3 # index = 4 is the method, index = 3 is the tpl, 
        subset_nmax = 12 # max num of subset shown in figure
        data_label_subset = [a[subset_index] for a in dataCalc]
        s = set(data_label_subset)
        data_subsets = []
        for a in s:
            #print (a)
            data_subsets.append([a, data_label_subset.count(a)])
        #sorted_hist = sorted(hist, key=lambda hist: hist[0])
        print("subsets = {0}".format(data_subsets))
    
        ######### plot
        xvect = range(0,n)
        yvect = range(0,n)  
        for k in range(0,n):
        	xvect[k] = dataRef[k][iRefZ]
    	yvect[k] = (dataCalc[k][1]-dataRef[k][iRefZ])/(1+dataRef[k][iRefZ])
        fig = pp.figure('Amazed output')
        ax = fig.add_subplot(111)
        color_ = ['r', 'g', 'b', 'y','c', 'm', 'y', 'k', 'r', 'g', 'b', 'y']
        #n=max(len(s),5)
        #color_=cm.rainbow(np.linspace(0,1,n))
        mylegend = [a[0] for a in data_subsets]
        for k in range(min(len(data_subsets),subset_nmax)):
            inds = [i for i,x in enumerate(data_label_subset) if x==data_subsets[k][0]]
    	#print("indices found are {0}".format(inds))
            xv = [xvect[i] for i in inds ]
            yv = [yvect[i] for i in inds ]
            ax.plot(xv, yv, 'x', label=mylegend[k], color=color_[k])
    
        #pp.legend((mylegend),loc=4)
        #pp.legend((mylegend))
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        pp.legend(loc=9, bbox_to_anchor=(0.5, -0.1))
        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Put a legend to the bottom of the current axis
        #pp.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
         #     fancybox=True, shadow=True, ncol=5)
        #pp.hist(yvect, 5000, normed=1, histtype='step', cumulative=True)
        ##bar
        #ind = np.arange(len(OY))
        #pp.plot(xvect, yvect, 'x')
        #ax.set_xscale('log')
        pp.grid(True) # Affiche la grille
        #pp.legend(('cos','sin'), 'upper right', shadow = True)
        pp.ylabel('(zcalc-zref)/(1+zref)')
        pp.xlabel('z reference')
        #pp.title('Amazed performance stats') # Titre
        pp.savefig( os.path.dirname(os.path.abspath(outFile)) + '/' +'diff.png', bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
        #pp.show()
        zoomedRange = 0.005
        ax.set_ylim([-zoomedRange, zoomedRange])
        #pp.show()
        pp.savefig( os.path.dirname(os.path.abspath(outFile)) + '/' +'diff_yzoomed.png', bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png


def loadRef(fname, reftype):
    """
    load the ref file data
    """ 
    print("loadref with type={}".format(reftype))
    dataArray = []
    f = open(fname)
    for line in f:
        lineStr = line.strip()
        if not lineStr.startswith('#'):
            #print lineStr
            data = lineStr.split("\t")
            if len(data)<2:
                data = lineStr.split(" ")
            data = [r for r in data if r != '']
            #print len(data)
            if reftype=="pfs":
                if(len(data) >= 7): #PFS
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = str(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = float(data[6])
                    d = [d0, d1, d2, d3, d4, d5, d6]
                    dataArray.append(d) 
            elif reftype=="muse":
                if(len(data) >= 2): #MUSE
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d = [d0, d1]
                    dataArray.append(d) 
            elif reftype=="keck":
                if(len(data) >= 7): #KECK
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = str(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = float(data[6])
                    d = [d0, d1, d2, d3, d4, d5, d6]
                    dataArray.append(d) 
            elif reftype=="simulm":
                if(len(data) >= 8): #SIMULM
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = str(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = float(data[6])
                    d7 = float(data[7])
                    d = [d0, d1, d2, d3, d4, d5, d6, d7]
                    dataArray.append(d) 
            elif reftype=="simueuclid2016":
                if(len(data) >= 8): #SIMUEuclid2016
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = str(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = float(data[6])
                    d7 = float(data[7])
                    d = [d0, d1, d2, d3, d4, d5, d6, d7]
                    dataArray.append(d) 
            elif reftype=="vvds2":
                if(len(data) >= 7): #VVDS2
                    d0 = str(data[0])
                    d1 = str(data[1])
                    d2 = float(data[2])
                    d3 = float(data[3])
                    d4 = float(data[4])
                    d5 = float(data[5])
                    d6 = float(data[6])
                    d = [d0, d1, d2, d3, d4, d5, d6]
                    dataArray.append(d) 
            elif reftype=="vuds":
                if(len(data) >= 6):
                    d0 = str(data[0])
                    d1 = float(data[1])
                    d2 = float(data[2])
                    d3 = float(data[3])
                    d4 = int(data[4])
                    d5 = float(data[5])
                    d = [d0, d1, d2, d3, d4, d5]
                    dataArray.append(d) 
            
    f.close()
    return dataArray

def loadCalc(fname):
    """
    load the calc redshift.csv data from file
    1 - Spectrum name
    2 - redshift
    3 - Merit
    4 - Method
    """ 
    dataArray = []
    f = open(fname)
    for line in f:
        lineStr = line.strip()
        if not lineStr.startswith('#'):
            #print lineStr
            data = lineStr.split("\t")
            data = [r for r in data if r != '']
            #print len(data)
            if(len(data) == 4): #linematching for example
                d0 = str(data[0])
                d1 = float(data[1])
                d2 = float(data[2])
                d3 = str(data[3])
                d = [d0, d1, d2, d3]
                dataArray.append(d) 
            if(len(data) >= 5): #chisquarer for example
                d0 = str(data[0])
                d1 = float(data[1])
                d2 = float(data[2])
                d3 = str(data[3])
                d4 = str(data[4])
                d = [d0, d1, d2, d3, d4]
                dataArray.append(d) 
    f.close()
    return dataArray
    
def loadDiff(fname):
    """
    load the diff diff.txt data from file

    """ 
    dataArray = []
    f = open(fname)
    for line in f:
        lineStr = line.strip()
        if not lineStr.startswith('#'):
            #print lineStr
            data = lineStr.split("\t")
            data = [r for r in data if r != '']
            #print len(data)
            if(len(data) == 15): #Spectrum ID	MAGI	ZREF	ZFLAG	ZCALC	MERIT	TPL	METHOD    SNR	SFR E(B-V) Sigma LogHalpha ExtValue DIFF
                d0 = str(data[0])
                d1 = float(data[1])
                d2 = float(data[2])
                d3 = float(data[3])
                d4 = float(data[4])
                d5 = float(data[5])
                d6 = str(data[6])
                d7 = str(data[7])
                d8 = float(data[8])
                d9 = float(data[9])
                d10 = float(data[10])
                d11 = float(data[11])
                d12 = float(data[12])
                d13 = float(data[13])
                d14 = float(data[14])
                d = [d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14]
                dataArray.append(d) 
    f.close()
    return dataArray

def ProcessStats( fname, zRange, magRange,  sfrRange, enablePlot = False ):
    data = loadDiff( fname );
    
    
    _baseOutputDirectory = os.path.dirname(os.path.abspath(fname))
    outputDirectory = os.path.join(_baseOutputDirectory, "stats_magmin{}magmax{}_zmin{}zmax{}_sfrmin{}sfrmax{}".format(magRange[0], magRange[1], zRange[0], zRange[1], sfrRange[0], sfrRange[1]))
    if os.path.exists( outputDirectory ) == False :        
        print("makedir: Output dir: "+outputDirectory)
        os.mkdir( outputDirectory )        
        
    n = (len(data))
    n2 = len(data[0])
    print "INFO: processing stats: n=" + str(n) + ", n2=" + str(n2)
    
    zmin = zRange[0]
    zmax = zRange[1]
    magmin = magRange[0]
    magmax = magRange[1]
    sfrmin = sfrRange[0]
    sfrmax = sfrRange[1]
    print("INFO: processing stats: zmin = {}, zmax = {}".format(zmin, zmax))
    print("INFO: processing stats: magmin = {}, magmax = {}".format(magmin, magmax))
    print("INFO: processing stats: sfrmin = {}, sfrmax = {}".format(sfrmin, sfrmax))
        
    indsForHist = []
    for x in range(0,n):
        if 0: #filter by flag
            flag = data[x][3]
            if flag!=2.0:
                continue
        if 0: #filter by n valid lines
            extvect = (data[x][13])
            if extvect != 2:
                continue
        zreference = (data[x][2])
        mag = (data[x][1])
        sfr = (data[x][9])
        if zreference <= zmax and zreference >= zmin and mag >= magmin and mag <= magmax and sfr >= sfrmin and sfr <= sfrmax:
            indsForHist.append(x)
            
    nSelected = len(indsForHist)
    namevect = range(0,nSelected)
    yvect = range(0,nSelected)
    mvect = range(0,nSelected)
    snrvect = range(0,nSelected)
    extvect = range(0,nSelected)
    sfrvect = range(0,nSelected)
    ebmvvect = range(0,nSelected)
    sigmavect = range(0,nSelected)
    fhalphavect = range(0,nSelected)
    zref = range(0,nSelected)
    zcalc = range(0,nSelected)
    for x in range(0,nSelected):
        namevect[x] = data[indsForHist[x]][0]
        
        mvect[x] = (data[indsForHist[x]][1])
        snrvect[x] = (data[indsForHist[x]][8])
        extvect[x] = (data[indsForHist[x]][13])
        sfrvect[x] = (data[indsForHist[x]][9])
        ebmvvect[x] = (data[indsForHist[x]][10]) 
        sigmavect[x] = (data[indsForHist[x]][11])
        fhalphavect[x] = (data[indsForHist[x]][12]) 
        zref[x] = (data[indsForHist[x]][2])
        zcalc[x] = (data[indsForHist[x]][4])
        
        yvect[x] = abs(data[indsForHist[x]][n2-1])/(1+zref[x])
        
    # export the filtered list
    foutpath = outputDirectory + '/' + 'stats_subset_list.txt'
    fout = open( foutpath, "w" )  
    outStr = ""
    outStr = outStr + "ID" + "\t"
    outStr = outStr + "ZCALC" + "\t"
    outStr = outStr + "DIFF" + "\t"
    outStr = outStr + "MAG" + "\t"
    outStr = outStr + "SNR" + "\t"
    outStr = outStr + "SFR" + "\t"
    outStr = outStr + "EBmV" + "\t"
    outStr = outStr + "Sigma" + "\t"
    outStr = outStr + "ZREF" + "\t"
    fout.write( outStr  + '\n')
    for x in range(0,nSelected):
        outStr = ""
        outStr = outStr + str(namevect[x]) + "\t"
        outStr = outStr + str(zcalc[x]) + "\t"
        outStr = outStr + str(yvect[x]) + "\t"
        outStr = outStr + str(mvect[x]) + "\t"
        outStr = outStr + str(snrvect[x]) + "\t"
        outStr = outStr + str(sfrvect[x]) + "\t"
        outStr = outStr + str(ebmvvect[x]) + "\t"
        outStr = outStr + str(sigmavect[x]) + "\t"
        outStr = outStr + str(zref[x]) + "\t"
        fout.write( outStr  + '\n')
    fout.close()
    # ...
    
    # plot the filtered list scatter zrelerr = f(1+zref)
    relerrorvect = range(0,nSelected)
    for x in range(0,nSelected):
        relerrorvect[x] = (zcalc[x] - zref[x])/(1+zref[x])
    fig = pp.figure('filtered dataset z relative error')
    pp.plot(zref, relerrorvect, 'x')
    pp.grid(True) # Affiche la grille
    #pp.legend(('cos','sin'), 'upper right', shadow = True)
    pp.ylabel('(zcalc-zref)/(1+zref)')
    pp.xlabel('z reference')
    pp.title('All spectra included') # Titre
    pp.savefig( os.path.join(outputDirectory, 'filteredset_relzerr.png'), bbox_inches='tight')
    # plot zoom 1
    zoomedRange = 0.005
    nmissing = 0
    for x in range(0,nSelected):
        if abs(relerrorvect[x])>zoomedRange:
            nmissing += 1
    pp.ylim([-zoomedRange, zoomedRange])
    #pp.show()
    spcStr = "spectrum"
    if nmissing>1:
        spcStr = "spectra"
    pp.title('{}/{} {} outside displayed range'.format(nmissing,nSelected, spcStr)) # Titre
    pp.savefig( os.path.join(outputDirectory, 'filteredset_relzerr_zoom_1.png'), bbox_inches='tight')
    # plot zoom 2
    zoomedRange = 0.0005
    nmissing = 0
    for x in range(0,nSelected):
        if abs(relerrorvect[x])>zoomedRange:
            nmissing += 1
    pp.ylim([-zoomedRange, zoomedRange])
    #pp.show()
    spcStr = "spectrum"
    if nmissing>1:
        spcStr = "spectra"
    pp.title('{}/{} {} outside displayed range'.format(nmissing,nSelected,spcStr)) # Titre
    pp.savefig( os.path.join(outputDirectory, 'filteredset_relzerr_zoom_2.png'), bbox_inches='tight')
    if enablePlot:
        pp.show()
        
    # plot the filtered list scatter zcalc = f(Zref)    
    fig = pp.figure('filtered dataset zcalc=f(zref)')
    pp.plot(zref, zcalc, 'x')
    pp.grid(True) # Affiche la grille
    #pp.legend(('cos','sin'), 'upper right', shadow = True)
    pp.ylabel('zcalc')
    pp.xlabel('zref')
    pp.title('All spectra included') # Titre
    pp.savefig( os.path.join(outputDirectory, 'filteredset_zcalc-vs-zref.png'), bbox_inches='tight')
        
    # plot Histograms
    # ******* large bins histogram
    vectErrorBins = [0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 1.0, 10.0]
    print 'the rough bins are: ' + str(vectErrorBins)
    foutpath = outputDirectory + '/' + 'stats_brief.txt'
    processHistogram( yvect, vectErrorBins, foutpath)

    # ******* fine bins histogram
    vectErrorBins = np.logspace(-5, 1, 50, endpoint=True)
    print 'the fine bins are: ' + str(vectErrorBins)
    foutpath = outputDirectory + '/' + 'stats.txt'
    processHistogram( yvect, vectErrorBins, foutpath)

    # ******* plot hist
    outFigFile = outputDirectory + '/' +'stats_hist.png'
    plotHist(yvect, outFigFile)

    nPercentileDepth = 2
    # ******* plot mag hist
    if 1:
        print("Plotting versus Mag")
        outFileNoExt = 'stats_versusMag_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, mvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='MAG', nPercentileDepth=nPercentileDepth)

    # ******* plot snr hist       
    if 1:
        print("Plotting versus Noise")
        outFileNoExt = 'stats_versusNoise_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, snrvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='SNR', nPercentileDepth=nPercentileDepth) 

    # ******* plot redshift hist       
    if 1:
        print("Plotting versus redshift")
        outFileNoExt = 'stats_versusRedshift_hist'
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt) 
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, zref, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='REDSHIFT', nPercentileDepth=nPercentileDepth) 
        
    # ******* plot sfr hist       
    if 1:
        print("Plotting versus sfr")
        outFileNoExt = 'stats_versusSFR_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, sfrvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='SFR', nPercentileDepth=nPercentileDepth) 

    # ******* plot EBmV hist       
    if 1:
        print("Plotting versus ebmv")
        outFileNoExt = 'stats_versusEBMV_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, ebmvvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='EBMV', nPercentileDepth=nPercentileDepth) 

    # ******* plot Sigma hist       
    if 1:
        print("Plotting versus sigma")
        outFileNoExt = 'stats_versusSigma_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, sigmavect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='SIGMA', nPercentileDepth=nPercentileDepth) 

    # ******* plot logHalpha hist       
    if 1:
        print("Plotting versus logfhalpha")
        outFileNoExt = 'stats_versusLogFHAlpha_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, fhalphavect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='LOGFHALPHA', nPercentileDepth=nPercentileDepth) 


    # ******* plot ext hist       
    if 1:
        print("Plotting versus ExtValue")
        outFileNoExt = 'stats_versusExt_hist' 
        outFilepathNoExt = os.path.join(outputDirectory,outFileNoExt)
        outdir = outputDirectory
        lstats.PlotAmazedVersusBinsHistogram(yvect, extvect, outdir, outFilepathNoExt, enablePlot=enablePlot, enableExport=1, mtype='EXT', nPercentileDepth=nPercentileDepth) 


    print '\n'

def processPerformance( fname, opt_preset ):
    dataDiff = loadDiff( fname );
    
    _baseOutputDirectory = os.path.dirname(os.path.abspath(fname))
    outputDirectory = os.path.join(_baseOutputDirectory, "performances".format())
    if os.path.exists( outputDirectory ) == False :        
        print("makedir: Output dir: "+outputDirectory)
        os.mkdir( outputDirectory )
        
    lperf.exportPerformances(dataDiff, outputDirectory, preset=opt_preset)

def processHistogram(yvect, bins, outFile=""):
    n = len(yvect)
    nbins = len(bins)
    #print 'number of bins = ' + str(nbins) + '\n'
    ybins = np.zeros(nbins)

    #print range(nbins-2, -1, -1)
    for x in range(0,n):
	for ex in range(0, nbins):
	#for ex in range(nbins-2, -1, -1):
	    if yvect[x]<bins[ex]:
	   	ybins[ex] = ybins[ex]+1
		#break;
    for ex in range(0, nbins):
        ybins[ex] = ybins[ex]/float(n)*100.0
        
    if not outFile=="":
        # convert to percentage, and save in file 
        fout = open( outFile, "w" )  
        for ex in range(0, nbins):
            outStr = str(ex) + '\t' + str(bins[ex]) + '\t' +str(ybins[ex])
            print outStr
            fout.write( outStr  + '\n')
        fout.close()
    return ybins

def plotHist(yvect, outFigFile):
    subsets_enable = False
    ######### plot
    fig = pp.figure('Amazed performance stats')
    ax = fig.add_subplot(111)
    #ax.plot(vectErrorBins, yVectErrorBins)
    if 0:
        pp.hist(yvect, 5000, range = (1e-5, 10), normed=1, histtype='step', cumulative=True)
    else:
        vectErrorBins = np.logspace(-5, 1, 500, endpoint=True)
        ybins = processHistogram(yvect=yvect, bins=vectErrorBins, outFile="")
        pp.plot(vectErrorBins, ybins, "+-")
        pp.xlim([1e-5, 10])
        pp.ylim([0, 100])
    ##bar
    #ind = np.arange(len(OY))
    #pp.plot(xvect, yvect, 'x')
    ax.set_xscale('log')
    pp.grid(True) # Affiche la grille
    #pp.legend(('cos','sin'), 'upper right', shadow = True)
    ylabel = 'Success Rate (percentage) (over {} spectra)'.format(len(yvect))
    pp.ylabel(ylabel)
    pp.xlabel('abs( (zcalc-zref)/(1+zref) )')
    #pp.title('Amazed performance stats') # Titre
    pp.savefig( outFigFile, bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
    #pp.show()
    
    if subsets_enable == False:    
        return        

    #*************** create subsets
    subset_index = 7 # index = 7 is the method
    subset_nmax = 8 # max num of subset shown in figure
    data_label_subset = [a[subset_index] for a in data]
    s = set(data_label_subset)
    data_subsets = []
    for a in s:
        #print (a)
        data_subsets.append([a, data_label_subset.count(a)])
    #sorted_hist = sorted(hist, key=lambda hist: hist[0])
    print("subsets = {0}".format(data_subsets))

    ######### plot the subsets
    fig = pp.figure('Amazed performance stats by subsets')
    ax = fig.add_subplot(111)
    mylegend = [a[0]+" (n= " + str(a[1])+")" for a in data_subsets]
    color_ = ['r', 'g', 'b', 'y','c', 'm', 'y', 'k']
    for k in range(min(len(data_subsets),subset_nmax)):
        inds = [i for i,x in enumerate(data_label_subset) if x==data_subsets[k][0]] 
        yv = [yvect[i] for i in inds ]
        pp.hist(yv, 5000, normed=1, histtype='step', cumulative=True, label=mylegend[k], color=color_[k])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    pp.legend(loc=9, bbox_to_anchor=(0.5, -0.1))
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))    
    ax.set_xscale('log')
    pp.grid(True) # Affiche la grille
    #pp.legend(('cos','sin'), 'upper right', shadow = True)
    pp.ylabel('Cumulative Histogram')
    pp.xlabel('abs(diff)')
    #pp.title('Amazed performance stats') # Titre
    pp.savefig( os.path.dirname(os.path.abspath(fname)) + '/' +'stats_hist_subsets.png', bbox_inches='tight') # sauvegarde du fichier ExempleTrace.png
    #pp.show()
    

def ProcessFailures( fname, fnameFailures ): 
    fFailures = open(fnameFailures, 'w')

    
    if 0: 
        f = open(fname, 'r')
        dataStr = f.read()
        f.close()
        data= ascii.read(dataStr, data_start=0, delimiter="\t")
        print data[0]
    else:
        data = loadDiff( fname );

    n = (len(data))
    n2 = len(data[0])
    print "INFO: processing failures: n=" + str(n) + ", n2=" + str(n2)
    
    f = open(fname, 'r')
    dataHeader = f.readline()
    f.close()
    if "#" in dataHeader:
        fFailures.write( dataHeader )
    #print range(nbins-2, -1, -1)
    for x in range(0,n):
        if abs(data[x][n2-1]) > 1e-2:
            for x2 in range(0,n2):
                fFailures.write( str(data[x][x2]) + "\t" )
            fFailures.write( "\n" )


    fFailures.close()
    print '\n'    

def ProcessFailuresSeqFile( fname, refFile, fnameFailuresSeqFile, fnameFailureRefFile ): 
    """
    1. Creates a seq file "fnameFailuresSeqFile" that can be used as an input .spectrumlist file for Amazed.
    2. Also creates the appropriate ref file in order to compare the results if a reprocess on this sub-list is carried out...
    """
    fFailures = open(fnameFailuresSeqFile, 'w')
    if 0: 
        f = open(fname, 'r')
        dataStr = f.read()
        f.close()
        data= ascii.read(dataStr, data_start=0, delimiter="\t")
        print data[0]
    else:
        data = loadDiff( fname );

    n = (len(data))
    n2 = len(data[0])
    print "INFO: processing failures Seq file: n=" + str(n) + ", n2=" + str(n2)
    
    f = open(fname, 'r')
    dataHeader = f.readline()
    f.close()
    if "#" in dataHeader:
        #fFailures.write( dataHeader )
        pass
    #print range(nbins-2, -1, -1)
    dataFailures_names_raw = []
    for x in range(0,n):
        if abs(data[x][n2-1]) > 1e-2:
            indName = 0;    
            spcName = str(data[x][indName]) + ".fits"
            dataFailures_names_raw.append(spcName)
            spcStr = "_atm_clean"
            noiseStr = "_noise"
            if spcName.find(spcStr)==-1:
                spcStr = "_F_"
                noiseStr = "_ErrF"
            noiseName = spcName.replace(spcStr, noiseStr)
            fFailures.write( spcName + "\t" )
            fFailures.write( noiseName )
            fFailures.write( "\n" )
    fFailures.close()
    
    
    #*************** open ref file
    fref = open(refFile, 'r')
    dataRefStr = fref.read()
    fref.close()
    dataRef = ascii.read(dataRefStr)
    dataRef_names = [a[0] for a in dataRef] 
    
    nref2 = len(dataRef[0])
    #*************** reorder ref data by filename
    inds = []
    fFailuresRefFile = open(fnameFailureRefFile, 'w')
    #readonlyN = 300;
    for iref,s in enumerate(dataRef_names):#[0:readonlyN]:
        #remove extension .fits from pandora results
        #s = s[0:-5]
        p = [i for i,x in enumerate(dataFailures_names_raw) if str(s) in x]
        if len(p) >0:
            #print "OK : index found : ref={0} and calc={1}".format(s,p)
            inds.append(p[0])
            for x2 in range(0,nref2):
                fFailuresRefFile.write( str(dataRef[iref][x2]) + "\t" )
            fFailuresRefFile.write( "\n" )
        else:
            inds.append(-1)
            print "ERROR : index not found for failure-ref file creation: {0}".format(s)
            #stop  
    
    
    fFailuresRefFile.close()
    
    
    
    
    print '\n'

def StartFromCommandLine( argv ) :
    global subsets_enable	
    usage = """usage: %prog [options]
    ex: python ./processAmazedOutputStats.py --ref=referenceRedshifts.txt --calc=amazedRedshifts.txt --type='vvds2' """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-r", u"--ref", help="reference redshift values",  dest="refFile", default="referenceRedshifts.txt")
    parser.add_option(u"-c", u"--calc", help="calculated redshift values",  dest="calcFile", default="output.txt")
    parser.add_option(u"-t", u"--type", help="reference redshift values type, choose between 'vvds1', 'vvds2' or 'pfs', 'vuds', 'simulm', 'simueuclid2016'",  dest="type", default="vvds")
    parser.add_option(u"-d", u"--diff", help="diff file (overrides the use of --ref, --calc and --type)",  dest="diffFile", default="")
    
    parser.add_option(u"-m", u"--magRange", help="magnitude range filter for the histograms",  dest="magRange", default="0.0 50.0")
    parser.add_option(u"-s", u"--sfrRange", help="sfr range filter for the histograms",  dest="sfrRange", default="-1.0 10000.0")
    parser.add_option(u"-z", u"--zRange", help="redshift range filter for the histograms",  dest="zRange", default="-1.0 20.0")

    parser.add_option(u"-p", u"--perfpreset", help="performance matrix preset, choose between 'simulm201606', 'simueuclid2016'",  dest="perfpreset", default="simulm201606")
    
    
    parser.add_option(u"-l", u"--computeLvl", help="compute level, choose between 'brief' or 'full'",  dest="computeLevel", default="brief")
    (options, args) = parser.parse_args()

    print "\n"
    if os.path.isdir(options.calcFile):
        print("Error: invalid input calc. file (directory)")
        exit()
        
    if os.path.isdir(options.refFile):
        print("Error: invalid input ref. file (directory)")
        exit()
        


    subsets_enable = False;

    if( len( args ) == 0 ) :
        if options.diffFile == "":
            filenameDiff = "diff.txt"
            filenameFailures = "failures.txt"
            filenameFailuresSeqFile = "failures.spectrumlist"
            filenameFailuresRefFile = "failures_ref.txt"
            outputPath =  os.path.dirname(os.path.abspath(options.calcFile)) + '/' + "stats/"
            outputFullpathDiff = outputPath + filenameDiff
            outputFullpathFailures = outputPath + filenameFailures
            outputFullpathFailuresSeqFile = outputPath + filenameFailuresSeqFile
            outputFullpathFailuresRefFile = outputPath + filenameFailuresRefFile
            if os.path.isdir(outputPath)==False:
                os.mkdir( outputPath, 0755 );
    
            if options.type == 'vvds1':
    	       print "Info: Using VVDS1 reference data file type"
    	       setVVDSRefFileType()
            elif options.type == 'vvds2':
                print "Info: Using VVDS2 reference data file type"
                setVVDS2RefFileType()
            elif options.type == 'pfs':
                print "Info: Using PFS reference data file type"
                setPFSRefFileType()
            elif options.type == 'keck':
                print "Info: Using KECK reference data file type"
                setKeckRefFileType()
            elif options.type == 'muse':
                print "Info: Using MUSE reference data file type"            
                setMuseRefFileType()
            elif options.type == 'vuds':
                print "Info: Using VUDS reference data file type"            
                setVUDSRefFileType()
            elif options.type == 'simulm':
                print "Info: Using SimuLM reference data file type"            
                setSIMULMRefFileType()
            elif options.type == 'simueuclid2016':
                print "Info: Using SimuEuclid2016 reference data file type"            
                setSIMUEuclid2016RefFileType()
            else:
                print("Info: No reference file type given (--type), using vvds by default.")
    
            ProcessDiff( options.refFile, options.calcFile, outputFullpathDiff, options.type )
            ProcessFailures( outputFullpathDiff, outputFullpathFailures)
            ProcessFailuresSeqFile( outputFullpathDiff, options.refFile, outputFullpathFailuresSeqFile, outputFullpathFailuresRefFile)
        else:
            if os.path.exists(options.diffFile):
                outputFullpathDiff = options.diffFile
            else:
                print("ERROR: Cannot find diff file, aborting: {}".format(options.diffFile))
        
        if  options.computeLevel == "full" or options.computeLevel == "hist":
            zRange = [-1.0, 20.0]
            zRange[0] = float(options.zRange.split(" ")[0])
            zRange[1] = float(options.zRange.split(" ")[1])
            magRange = [0.0, 40.0]
            magRange[0] = float(options.magRange.split(" ")[0])
            magRange[1] = float(options.magRange.split(" ")[1])
            sfrRange = [0.0, 10000.0]
            sfrRange[0] = float(options.sfrRange.split(" ")[0])
            sfrRange[1] = float(options.sfrRange.split(" ")[1])        
            ProcessStats( outputFullpathDiff, zRange, magRange, sfrRange )
        
        if  options.computeLevel == "full" or options.computeLevel == "perf":        
            processPerformance( outputFullpathDiff, opt_preset = options.perfpreset )
            
    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()
    
if __name__ == '__main__':
    Main( sys.argv )

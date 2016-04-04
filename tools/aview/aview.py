# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 11:26:41 2015

@author: aschmitt
"""

__all__ = ['plotRes']

import sys
import os
import optparse

import resparser as rp
import aviewplot as avp
import resultstat as rstat
import chisquare as chisq

def plotRes(resDir, spcName, tplpath, redshift, iextremaredshift, diffthres, failureIdx, enablePlot=True):
    print('using amazed results full path: {0}'.format(resDir))
    s = rp.ResParser(resDir)
    print(s) 
    
    if spcName == "":
        print("Using Diffthreshold: {}".format(diffthres))
        resList = rstat.ResultList(resDir, diffthreshold=float(diffthres), opt='brief')        
        print("Results, N found failures = {0}".format(resList.n))
        indice = int(failureIdx)
        if indice < resList.n:
            spcName = resList.list[indice].name
            print("\n")
            print("plotRes: spc = {}".format(spcName))
            print("plotRes:, zref = {}".format(resList.list[indice].zref))
            print("plotRes:, zcalc = {}".format(resList.list[indice].zcalc))
            print("plotRes:, diff = {}".format(resList.list[indice].zdiff))
            print("\n")

        else:
            print("WARNING: index is not within the failure list size range... aborting".format())
            return
        
    spath = s.getSpcFullPath(spcName)
    print('Spc path is: {0}'.format(spath))

    #prepare exportPath if necessary
    if enablePlot==False:
        displaySpcPath = os.path.join(s.getDisplaysDirPath(), spcName)
        if not os.path.exists(displaySpcPath):
            os.makedirs(displaySpcPath)
    else:
        displaySpcPath = ""

    # PLOT Chisquare 
    spctag = spcName;
    extensions = ['.fits', '.FITS', '.txt', '.TXT', '.dat', '.DAT']
    for e in extensions:
        if spctag.endswith(e):
            iend = len(e)
            spctag = spctag[:-iend]
            break
    path = os.path.join(resDir, spctag)
    
    [chipathlist, chinamelist] = s.getAutoChi2FullPath(spcName)
    print("chipathlist found : {}".format(chipathlist))       
       
    for ii,chipath in enumerate(chipathlist):
        print('Trying to plot chi2 using full path: {}'.format(chipath))
        if os.path.exists(chipath):
            print('plot chi2 using full path: {}'.format(chipath))
            chi = chisq.ResultChisquare(chipath, stype=os.path.splitext(chinamelist[ii])[0])
            #print(chi) 
            chi.plot(showContinuumEstimate=False, showExtrema=True, showAmbiguities=False, enablePlot=enablePlot, exportPath=displaySpcPath)


    try:
        if not iextremaredshift == "": 
            idxExtrema = int(iextremaredshift)
        else:
            idxExtrema = int(0)
        print('using idxExtrema: {}'.format(idxExtrema))
            

        tplpath = s.getAutoTplFullPath(spcName, idxExtrema)
        print("full extrema model path: {}".format(tplpath))
        if not os.path.exists(tplpath):
            print("reset tplpath = {}".format(""))
            tplpath = ""
            forceTplAmplitude = -1
            forceTplDoNotRedShift = 0
        else:
            forceTplAmplitude = 1
            forceTplDoNotRedShift = 1
    except:
        pass
            

    npath = ""
#    spctag = "atm_clean" #for VVDS
#    if spcName.find(spctag)!=-1:
#        noiseName = spcName.replace("atm_clean", "noise")
#        npath = s.getNoiseFullPath(noiseName)
#    spctag = "-W-F_" #for PFS
#    if spcName.find(spctag)!=-1:
#        noiseName = spcName.replace(spctag, "-W-ErrF_")
#        npath = s.getSpcFullPath(noiseName)
#        
    npath = s.getNoiseFullPath(spcName)
    if not os.path.exists(npath):
        npath = ""
    print('Noise path is: {}'.format(npath))

    
    if redshift == "":
        if not iextremaredshift == "": 
            idxExtrema = int(iextremaredshift)
        else:
            idxExtrema = int(0)
        print('using idxExtrema: {}'.format(idxExtrema)) 
        if 'chi' in locals():
            zval = chi.getExtrema(idxExtrema)
        
        if not 'zval' in locals():
            print("WARNING: could not use first extrema as redshift value !!!!")
            zval = s.getRedshiftVal(spcName)
        #overrride zval
        #zval = 1.1164
    else:
        print("WARNING: using custom user defined redshift value !!!!")
        zval = float(redshift)    
    
    print('Redshift is: {0}'.format(zval))
    #print(s.getRedshiftTpl(spcName))
    
    if tplpath == "":
        tpath = s.getTplFullPath(s.getRedshiftTpl(spcName)) #AUTO from amazed results
        try:
            forceTplAmplitude = chi.amazed_fitamplitude[idxExtrema]
        except:
            forceTplAmplitude = -1
        
        if 0:
            #tpath = ""
            #tpath = s.getTplFullPath("NEW_Im_extended.dat") #itpl = 0
            #tpath = s.getTplFullPath("NEW_Im_extended_blue.dat") #itpl = 1
            #tpath = s.getTplFullPath("NEW_Sbc_extended.dat") #itpl = 2
            #tpath = s.getTplFullPath("Scd.txt") #itpl = 3
            #tpath = s.getTplFullPath("StarBurst1.txt") #itpl = 4
            #tpath = s.getTplFullPath("StarBurst2.txt") #itpl = 5
            #tpath = s.getTplFullPath("StarBurst3.txt") #itpl = 6
            #tpath = s.getTplFullPath("BulgedataExtensionData.dat") #itpl = 7
            #tpath = s.getTplFullPath("EW_SB2extended.dat") #itpl = 8
            #tpath = s.getTplFullPath("EdataExtensionData.dat") #itpl = 10
            #tpath = s.getTplFullPath("NEW_E_extendeddataExtensionData.dat") #itpl = 10, 
            tpath = s.getTplFullPath("EllipticaldataExtensionData.dat") #itpl = 11
            #tpath = s.getTplFullPath("NEW_E_extendeddataExtensionData.dat") #itpl = 12
            #tpath = s.getTplFullPath("s0dataExtensionData.dat") #itpl = 13  
            #tpath = s.getTplFullPath("sadataExtensionData.dat")
            #tpath = s.getTplFullPath("vvds_reddestdataExtensionData.dat")
            forceTplAmplitude = -1
            forceTplAmplitude = 7e-5
        
        print('using forceTplAmplitude: {}'.format(forceTplAmplitude))
        
        forceTplDoNotRedShift = 0
    elif tplpath==-1:
        tpath = ""
    else:
        tpath = tplpath
    
    print('Tpl path is: {0}'.format(tpath) )  
    cpath = s.getCatalogFullPath()
    if 0:
        path = "/home/aschmitt/data/pfs/pfs_lbg/amazed/RayCatalogs"
        name = "raycatalogamazedvacuum.txt"
        cpath = os.path.join(path,name)
    
    print('Catalog path is: {0}'.format(cpath) )  
    
    # PLOT Spectrum/tpl/lines view
    avp1 = avp.AViewPlot(spath, npath, tpath, cpath, zval, forceTplAmplitude=forceTplAmplitude, forceTplDoNotRedShift = forceTplDoNotRedShift, enablePlot=enablePlot)
    if enablePlot==False:
        avp1.exportDisplays(displaySpcPath);

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./aview.py """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-x", u"--xaxis", help="force the use of indexes for the xaxis (--xaxis=""index""), or enable the use of calibrated wavelength (--xaxis=""angstrom"")",  dest="xAxisType", default="angstrom")
    parser.add_option(u"-d", u"--dir", help="path to the amazed results directory (/output/)",  dest="resDir", default="./output")
    parser.add_option(u"-s", u"--spc", help="name of the spectrum to be plotted",  dest="spcName", default="")
    parser.add_option(u"-t", u"--tpl", help="path of the template to be plotted",  dest="tplpath", default="")
    parser.add_option(u"-z", u"--redshift", help="z to be plotted",  dest="redshift", default="")
    parser.add_option(u"-e", u"--iextremaredshift", help="extrema index for the z to be plotted",  dest="iextremaredshift", default="")
    parser.add_option(u"-f", u"--failurediffthres", help="diff threshold for the extraction of the failure spectra",  dest="diffthres", default=-1)
    parser.add_option(u"-i", u"--failureindex", help="failure index",  dest="failureindex", default=0)
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        plotRes(options.resDir, options.spcName, options.tplpath, options.redshift, options.iextremaredshift, options.diffthres, options.failureindex)
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

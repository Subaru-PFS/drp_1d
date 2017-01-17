# -*- coding: utf-8 -*-
"""
"""
 
__all__ = ['export']

import sys
import os
import argparse
import re

import spectrum
import reference

def export(spcListPath, suffix, enableOverlay, othersuffix, otherpath, refpath):
    
    print('using spclist full path: {}'.format(spcListPath))
    f = open(spcListPath)
    spcList = [] 
    for line in f:
        lineStr = line.strip()
        if not lineStr.startswith('#'):
            #print lineStr
            #data = lineStr.split("\t| ")
            data = re.split("\t| ", lineStr)
            data = [r for r in data if r != '']
            if(len(data) ==2):
                _name  = data[0]
                _noise =  data[1]
                spcList.append([_name, _noise])
    print("spcList size = {}".format(len(spcList)))
    if len(spcList)<1:
         return
         
    #load catalog information if any
    enableRefFile = False
    if not refpath=="":
        if os.path.exists(refpath):
            refCatalog = reference.Reference(referencepath=refpath, rtype="pfs")
            enableRefFile = True

    basePath = os.path.split(spcListPath)[0]
    displaySpcPath = os.path.join(basePath, "display")
    if not os.path.exists(displaySpcPath):
        os.makedirs(displaySpcPath)
    
    for e in spcList:
        try:
            s1 = spectrum.Spectrum(e[0])
            if enableOverlay:
                eTFTmp = e[0].replace(suffix, othersuffix)
                eTFName = os.path.split(eTFTmp)[1]
                eTF = os.path.join(otherpath, eTFName)
                s2 = spectrum.Spectrum(eTF)
                spcExportPathFull = os.path.join(displaySpcPath, "{}.png".format(e[0]))
                print("exporting to : {}".format(spcExportPathFull))
                
                # some specifics for PFS sim2016
                #s1.applyWeight(1e-17)
                s2.applyLambdaCrop(3800, 12600)
                label1 = "Flux with 3H instrument signature"
                label2 = "True Flux"
                title_overide = ""
                if enableRefFile:
                    idxobj = refCatalog.findIdx(eTFName)
                    if not idxobj == -1:
                        title_override = refCatalog.getTag(idxobj)
                ######
        
                s1.plotCompare(s2, 1.0, modellinetype = "k-", exportPath=spcExportPathFull, label2=label2, label1=label1, title_suffix=title_override)
            else:
                s1.plot(saveFullDirPath = displaySpcPath, lstyle="b-")
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            print("Problem with spectrum")
            raise
            continue
                    
    f.close()
    
        

def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--spectrumlist", dest="spcList", default="./list.spectrumlist",
                    help="path to the amazed spectrumlist file")  
    parser.add_argument("-s", "--suffix", dest="suffix", default="_F.",
                    help="suffix of the main spectrum, as in the spectrumList file")
                    
    #overlay spectrum: _TF spectrum for example
    parser.add_argument("-a", "--addOverlay", dest="addOverlay", action='store_true',
                    help="enable overlay another spectrum to the existing one") 
    parser.add_argument("-o", "--otherSuffix", dest="othersuffix", default="_TF.",
                    help="suffix of the overlay spectrum, so that '_F.' will be replaced by this suffix")
    parser.add_argument("-p", "--otherPath", dest="otherpath", default="./",
                    help="path of the other spc to be overlaid") 
    
    #catalog information                
    parser.add_argument("-r", "--refPath", dest="refpath", default="",
                    help="path of the ref file for catalog information")    

    options = parser.parse_args()    
    print(options)
    enableOverlay = options.addOverlay
    
    if os.path.exists(options.spcList):
        export(options.spcList, options.suffix, enableOverlay, options.othersuffix, options.otherpath, refpath=options.refpath)
    else :
        print("Error: invalid input arguments")
        exit()


def Main( argv ) :	
    try:
        StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()
   
if __name__ == '__main__':
    Main( sys.argv )
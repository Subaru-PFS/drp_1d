# -*- coding: utf-8 -*-
"""
"""
 
__all__ = ['export']

import sys
import os
import argparse
import re

import spectrum

def export(spcListPath, suffix, othersuffix, otherpath):
    
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

    displaySpcPath = os.path.join("", "display")
    if not os.path.exists(displaySpcPath):
        os.makedirs(displaySpcPath)
    
    for e in spcList:
        try:
            s1 = spectrum.Spectrum(e[0])
            eTFTmp = e[0].replace(suffix, othersuffix)
            eTFName = os.path.split(eTFTmp)[1]
            eTF = os.path.join(otherpath, eTFName)
            s2 = spectrum.Spectrum(eTF)
            spcExportPathFull = os.path.join(displaySpcPath, "{}.png".format(e[0]))
            s1.plotCompare(s2, 1.0, modellinetype = "k-", exportPath=spcExportPathFull)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            continue
                    
    f.close()
    
        

def StartFromCommandLine( argv ) :	
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--spectrumlist", dest="spcList", default="./list.spectrumlist",
                    help="path to the amazed spectrumlist file")  
    parser.add_argument("-s", "--suffix", dest="suffix", default="_F.",
                    help="suffix of the main spectrum, as in the spectrumList file")
    parser.add_argument("-o", "--otherSuffix", dest="othersuffix", default="_TF.",
                    help="suffix of the overlay spectrum, so that '_F.' will be replaced by this suffix")
    parser.add_argument("-p", "--otherPath", dest="otherpath", default="./",
                    help="path of the other spc to be overlaid")    

    options = parser.parse_args()    
    print(options)
    
    if os.path.exists(options.spcList):
        export(options.spcList, options.suffix, options.othersuffix, options.otherpath)
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
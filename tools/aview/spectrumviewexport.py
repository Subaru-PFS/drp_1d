# -*- coding: utf-8 -*-
"""
"""
 
__all__ = ['export']

import sys
import os
import optparse
import re

import spectrum

def export(spcListPath, spcName=""):
    
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
        s1 = spectrum.Spectrum(e[0])
        eTF = e[0].replace("_F.", "_TF.")
        s2 = spectrum.Spectrum(eTF)
        spcExportPathFull = os.path.join(displaySpcPath, "{}.png".format(e[0]))
        s1.plotCompare(s2, 1.0, modellinetype = "k-", exportPath=spcExportPathFull)
        
                    
    f.close()
    
        

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./spectrumviewexport.py -d path_to_res_dir"""
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-i", u"--spectrumlist", help="path to the amazed spectrumlist file",  dest="spcList", default="./list.spectrumlist")
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        export(options.spcList)
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
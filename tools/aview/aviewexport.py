# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:20:56 2016

@author: aschmitt
"""
 
__all__ = ['export']

import sys
import os
import optparse

import resparser as rp
import resultstat as rstat
import aview

def export(resDir):
    
    print('using amazed results full path: {0}'.format(resDir))
    s = rp.ResParser(resDir)
    print(s) 
    
    diffthres = -1.0
    #print("Using Diffthreshold: {}".format(diffthres))
    resList = rstat.ResultList(resDir, diffthreshold=float(diffthres), opt='brief')        
    print("Results, N found failures = {0}".format(resList.n))
    for indice in range(resList.n):
        _spcName = ""
        _tplpath = ""
        _redshift = ""
        _iextremaredshift = 0
        _diffthres = diffthres
        _failureindex = indice
        _resDir = resDir
        aview.plotRes(_resDir, _spcName, _tplpath, _redshift, _iextremaredshift, _diffthres, _failureindex, enablePlot=False)
        

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./aviewexport.py -d path_to_res_dir"""
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-d", u"--dir", help="path to the amazed results directory (/output/)",  dest="resDir", default="./output")
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        export(options.resDir)
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
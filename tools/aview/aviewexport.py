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

def export(resDir, spcName, iextremaredshift):
    
    print('using amazed results full path: {0}'.format(resDir))
    s = rp.ResParser(resDir)
    print(s) 
    
    diffthres = -1.0
    #print("Using Diffthreshold: {}".format(diffthres))
    resList = rstat.ResultList(resDir, diffthreshold=float(diffthres), opt='brief')        
    print("Results, N found failures = {0}".format(resList.n))
    nDisp = resList.n
    if not spcName=="":
        nDisp = 1
    for indice in range(nDisp):
        _spcName = spcName
        _tplpath = ""
        _redshift = ""
        _iextremaredshift = iextremaredshift
        _diffthres = diffthres
        _failureindex = indice
        _resDir = resDir
        try:
            aview.plotRes(_resDir, _spcName, _tplpath, _redshift, _iextremaredshift, _diffthres, _failureindex, enablePlot=False)
        except:
            pass

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options]
    ex: python ./aviewexport.py -d path_to_res_dir"""
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-d", u"--dir", help="path to the amazed results directory (/output/)",  dest="resDir", default="./output")
    parser.add_option(u"-s", u"--spc", help="name of the spectrum to be plotted",  dest="spcName", default="")
    parser.add_option(u"-e", u"--iextremaredshift", help="extrema index for the z to be plotted",  dest="iextremaredshift", default="")
    (options, args) = parser.parse_args()

    if( len( args ) == 0 ) :
        export(options.resDir, options.spcName, options.iextremaredshift)
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
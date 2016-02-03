#!/bin/python
# -*- coding:utf-8 -*-

import os
from shutil import rmtree
import sys

amazedPath = os.path.expanduser ( "../../../amazed" )
amazedExecutable = amazedPath + "/bin/amazed-0.0.0"
cpfPath = os.path.expanduser ( "../.." )
successStrings = [ ]
failureStrings = [ ]
errorStrings = [ ]

def amazed ( amazedPathname, spectrumlist, templatedir, spectrumPath, linecatalogPathname, parametersPathname ):
    """
    Runs amazed.
    """
    amazedString = "{} --input {} --templatedir {} --spectrum {} -y {} -o output -m linematching2 --parameters {} --thread-count 1 > /dev/null".format (
        amazedPathname,
        spectrumlist,
        templatedir,
        spectrumPath,
        linecatalogPathname,
        parametersPathname )
    os.system ( amazedString )

def cleanup ( ):
    try:
        rmtree ( "output" )
    except:
        pass

def printErrors ( ):
    print ( "ErroredTests:" )
    for error in range ( len ( errorStrings ) ):
        print ( "{}".format ( errorStrings [ error ] ) )

def printFailures ( ):
    print ( "FailedTests:" )
    for failure in range ( len ( failureStrings ) ):
        print ( "{}".format ( failureStrings [ failure ] ) )

def printSuccesses ( ):
    print ( "SuccessfulTests:" )
    for success in range ( len ( successStrings ) ):
        print ( "{}".format ( successStrings [ success ] ) )

def test_good_input_good_output ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testGood.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/all_methods_dev.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 0.077164
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
        print ( "{}: redshift {} too different from 0.077164".format ( sys._getframe ( ).f_code.co_name, redshift ) )

def test_bad_input_bad_output ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testBad.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/all_methods_dev.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 3.59831
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )

if __name__ == "__main__":
    cleanup ( )
    try:
        test_good_input_good_output ( )
    except Exception as e:
        errorStrings.append ( "test_good_input_good_output: {}".format ( e ) )
    cleanup ( )
    try:
        test_bad_input_bad_output ( )
    except Exception as e:
        errorStrings.append ( "test_bad_input_bad_output: {}".format ( e ) )
    cleanup ( )

    printSuccesses ( )
    printFailures ( )
    printErrors ( )
    

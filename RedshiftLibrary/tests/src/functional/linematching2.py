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
    #amazedString = "{} --input {} --templatedir {} --spectrum {} -y {} -o output -m linematching2 --parameters {} --thread-count 1".format (
    amazedString = "{} --input {} --templatedir {} --spectrum {} -y {} -o output -m linematching2 --parameters {} --thread-count 1 > /dev/null 2>&1".format (
        amazedPathname,
        spectrumlist,
        templatedir,
        spectrumPath,
        linecatalogPathname,
        parametersPathname )
    try:
        return os.system ( amazedString )
    except Exception as e:
        return e

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

def wrapCall ( testFunctionName ):
    try:
        testFunctionName ( )
    except Exception as e:
        errorStrings.append ( "{}: {}".format ( testFunctionName.__name__, e ) )
    cleanup ( )
        
def testGoodSpectrumDynamicallyYieldsGoodRedshift ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testGood.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/linematching2Dynamic.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 0.077164
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
        print ( "{}: redshift {} too different from 0.077164".format ( sys._getframe ( ).f_code.co_name, redshift ) )

def testGoodSpectrumStaticallyYieldsGoodRedshift ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testGood.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/linematching2Static.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 0.077164
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
        print ( "{}: redshift {} too different from 0.077164".format ( sys._getframe ( ).f_code.co_name, redshift ) )

def testBadSpectrumStaticallyYieldsBadRedshift ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testBad.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/linematching2Static.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 3.59831
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )

def testBadSpectrumDynamicallyYieldsBadRedshift ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testBad.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/linematching2Dynamic.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 3.59831
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )

def testInvalidFluxFilenameOnSpectrumlistCausesAbortion ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testInvalid.spectrumlist"
    retval = amazed ( amazedExecutable,
                      spectrumlist,
                      cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
                      cpfPath + "/test/data/linematchingFunctional/batch6",
                      cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
                      cpfPath + "/test/data/linematchingFunctional/linematching2Static.json" )
    logFile = open ( "output/log.txt" )
    logString = logFile.read ( )
    failed = False
    if retval == 0:
        failed = True
        print ( "{}: Amazed returned a value of 0 to the operating system, despite input error.".format ( sys._getframe ( ).f_code.co_name ) )
    if not "Error: invalid file name 10000663000008vacLine_Invalid.fits specified in spectrumlist." in logString:
        failed = True
        print ( "{}: Amazed did not inform user of inexistent file listed on spectrumlist.".format ( sys._getframe ( ).f_code.co_name ) )

    if failed:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )

def testInvalidParameterNameOnJSONCausesAbortion ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testGood.spectrumlist"
    retval = amazed ( amazedExecutable,
                      spectrumlist,
                      cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
                      cpfPath + "/test/data/linematchingFunctional/batch6",
                      cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
                      cpfPath + "/test/data/linematchingFunctional/linematching2InvalidParameter.json" )
    logFile = open ( "output/log.txt" )
    logString = logFile.read ( )
    failed = False
    if retval == 0:
        failed = True
        print ( "{}: Amazed returned a value of 0 to the operating system, despite invalid parameter on JSON.".format ( sys._getframe ( ).f_code.co_name ) )
    if not "Error: invalid parameter specified in JSON." in logString:
        failed = True
        print ( "{}: Amazed did not inform user of invalid parameter present on JSON.".format ( sys._getframe ( ).f_code.co_name ) )

    if failed:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )

def testMissingParameterOnJSONCausesDefault ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testGood.spectrumlist"
    retval = amazed ( amazedExecutable,
                      spectrumlist,
                      cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
                      cpfPath + "/test/data/linematchingFunctional/batch6",
                      cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
                      cpfPath + "/test/data/linematchingFunctional/linematching2MissingParameter.json" )
    logFile = open ( "output/parameters.json" )
    logString = logFile.read ( )
    failed = False
    if retval != 0:
        failed = True
        print ( "{}: Amazed failed to reduce good data, when parameter missing on JSON.".format ( sys._getframe ( ).f_code.co_name ) )
    if not "\"minsize\": \"3\"" in logString:
        failed = True
        print ( "{}: Amazed did use default parameter value when parameter was missing in JSON.".format ( sys._getframe ( ).f_code.co_name ) )

    if failed:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )

def testInvalidParameterValueOnJSONCausesDefault ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testGood.spectrumlist"
    retval = amazed ( amazedExecutable,
                      spectrumlist,
                      cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
                      cpfPath + "/test/data/linematchingFunctional/batch6",
                      cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
                      cpfPath + "/test/data/linematchingFunctional/linematching2InvalidValue.json" )
    failed = False
    if retval != 0:
        failed = True
        print ( "{}: Amazed failed to reduce good data, when parameter has invalid value on JSON.".format ( sys._getframe ( ).f_code.co_name ) )
    try:
        logFile = open ( "output/parameters.json" )
        logString = logFile.read ( )
    except:
        failed = True
        print ( "{}: missing expected output/parameters.json file.".format ( sys._getframe ( ).f_code.co_name ) )
        logString = ""
    if "\"winsize\": \"-250\"" in logString:
        failed = True
        print ( "{}: Amazed did use invalid parameter value present in JSON.".format ( sys._getframe ( ).f_code.co_name ) )

    if failed:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )

def testLinetypeAbsorptionOnlyWorks ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testBad.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/linematching2AbsorptionOnly.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 1.46877
    if not expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
        print ( "{}: calculated redshift {} equal to emission-only redshift {}".format ( sys._getframe ( ).f_code.co_name, redshift, expectedRedshift ) )

def testLinetypeEmissionOnlyWorks ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testBad.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/linematching2EmissionOnly.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 1.46877
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
        print ( "{}: redshift {} too different from {}".format ( sys._getframe ( ).f_code.co_name, redshift, expectedRedshift ) )

def testLinetypeAllWorks ( ):
    spectrumlist = cpfPath + "/test/data/linematchingFunctional/batch6testBad.spectrumlist"
    amazed ( amazedExecutable,
             spectrumlist,
             cpfPath + "/test/data/linematchingFunctional/ExtendedGalaxyEL3",
             cpfPath + "/test/data/linematchingFunctional/batch6",
             cpfPath + "/test/data/linematchingFunctional/linecatalogamazedvacuum_B7C.txt",
             cpfPath + "/test/data/linematchingFunctional/linematching2AllLinetypes.json" )
    redshiftFile = open ( "output/redshift.csv" )
    redshiftString = redshiftFile.read ( )
    redshift = float ( redshiftString.split ( ) [ 1 ] )
    expectedRedshift = 1.46877
    if not expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        successStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
    else:
        failureStrings.append ( "{}".format ( sys._getframe ( ).f_code.co_name ) )
        print ( "{}: redshift {} too similar to emission only redshift {}".format ( sys._getframe ( ).f_code.co_name, redshift, expectedRedshift ) )

        
if __name__ == "__main__":
    cleanup ( )
    wrapCall ( testGoodSpectrumStaticallyYieldsGoodRedshift )
    wrapCall ( testGoodSpectrumDynamicallyYieldsGoodRedshift )
    wrapCall ( testBadSpectrumStaticallyYieldsBadRedshift )
    wrapCall ( testBadSpectrumDynamicallyYieldsBadRedshift )
    wrapCall ( testInvalidFluxFilenameOnSpectrumlistCausesAbortion )
    wrapCall ( testInvalidParameterNameOnJSONCausesAbortion )
    wrapCall ( testMissingParameterOnJSONCausesDefault )
    wrapCall ( testInvalidParameterValueOnJSONCausesDefault )
    wrapCall ( testLinetypeAbsorptionOnlyWorks )
    wrapCall ( testLinetypeEmissionOnlyWorks )
    wrapCall ( testLinetypeAllWorks )
            
    printSuccesses ( )
    printFailures ( )
    printErrors ( )

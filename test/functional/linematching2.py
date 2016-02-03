#!/bin/python
# -*- coding:utf-8 -*-

import os
from shutil import rmtree

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

def xmlHeader ( ):
    print ( "<?xml version=\"1.0\" encoding='ISO-8859-1' standalone='yes' ?>" )
    print ( "<TestRun>" )

def xmlFailures ( ):
    print ( "<FailedTests>" )
    for failure in range ( len ( failureStrings ) ):
        print ( "<Test id=\"{}\">".format ( failure + 1 ) )
        print ( "<Name>{}</Name>".format ( failureStrings [ failure ] ) )
        print ( "</Test>" )
    print ( "</FailedTests>" )

def xmlSuccesses ( ):
    print ( "<SuccessfulTests>" )
    for success in range ( len ( successStrings ) ):
        print ( "<Test id=\"{}\">".format ( success + 1 ) )
        print ( "<Name>{}</Name>".format ( successStrings [ success ] ) )
        print ( "</Test>" )
    print ( "</SuccessfulTests>" )

def xmlStatistics ( ):
    print ( "<Statistics>" )
    print ( "<Tests>{}</Tests>".format ( len ( successStrings ) + len ( failureStrings ) ) )
    print ( "<FailuresTotal>{}</FailuresTotal>".format ( len ( failureStrings ) ) )
    print ( "<Errors>{}</Errors>".format ( len ( errorStrings ) ) )
    print ( "<Failures>{}</Failures>".format ( len ( failureStrings ) ) )
    print ( "</Statistics>" )
    print ( "</TestRun>" )

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
    expectedRedshift = 3.59831
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        successStrings.append ( "test_good_input_good_output" )
    else:
        failureStrings.append ( "test_good_input_good_output" )

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
    expectedRedshift = 3.56089
    if expectedRedshift - 1e-4 <= redshift <= expectedRedshift + 1e-4:
        failureStrings.append ( "test_bad_input_bad_output" )
    else:
        successStrings.append ( "test_bad_input_bad_output" )

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

    xmlHeader ( )
    xmlSuccesses ( )
    xmlFailures ( )
    xmlStatistics ( )

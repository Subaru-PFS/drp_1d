#! /usr/bin/python

import urllib2
import sys
import optparse
import os
import tarfile
import shutil

thirdPartyDir = os.path.normpath(os.getcwd()+"/../thirdparty")+"/"

def ExtractTarGZ( tarPath, destPath ) :

    if os.path.exists( destPath ):
        print "File or folder: "+os.path.normpath( destPath )+" already exist, extraction skipped..."
        return False

    else :

        tfile = tarfile.open( tarPath, 'r:gz')

        print "Extracting: \n\tFrom: " + tarPath + " to: " + destPath

        extractDir = os.path.dirname(destPath)+"/"
        tfile.extractall( extractDir )
        extractedPath = os.path.normpath ( extractDir+tfile.getmembers()[0].name )

        os.rename( extractedPath, destPath )

        return True

def DownloadHTTPFile( fileUrl, localFilePath ) :

    if localFilePath == None :
        localFilePath = os.path.basename( fileUrl )

    # Check if file already exist
    if os.path.exists(localFilePath) :
        print "File: "+localFilePath+" already exist, download skipped..."
        return False

    # Do HTTP request
    try :
        urlfile = urllib2.urlopen( fileUrl )
    except urllib2.URLError as e:
        print "Download from: " + fileUrl + "failed.\nReason are:" + str( e.reason )
        return

    if urlfile.info().has_key('Content-Length'):
        totalFileSize = int (urlfile.info()['Content-Length'] )
    else :
        totalFileSize = 100

    localFile = open(localFilePath, 'wb')

    print "Downloading: \n\tFrom: " + fileUrl + "\n\tTo: " + localFilePath

    #Download file
    data_list = []
    chunk = 4096*10
    size = 0
    while 1:
        data = urlfile.read(chunk)
        if not data:
            print " OK !"
            break
        size += len( data )
        localFile.write(data)

        progress = float( size ) / totalFileSize * 100.0
        sys.stdout.write("\r[%f%%]" %progress )

        sys.stdout.flush()

    localFile.close()

    return True

libDict = {
    "cppunit": { "path": thirdPartyDir+"cppunit-1.12.1", "src": "http://downloads.sourceforge.net/cppunit/cppunit-1.12.1.tar.gz" },
    "boost": { "path":  thirdPartyDir+"boost-1.57.0", "src": "http://downloads.sourceforge.net/project/boost/boost/1.57.0/boost_1_57_0.tar.gz" },
    "gsl": { "path":  thirdPartyDir+"gsl-2.1", "src": "http://ftp.igh.cnrs.fr/pub/gnu/gsl/gsl-2.1.tar.gz" },
    "fftw": { "path":  thirdPartyDir+"fftw-3.3.3", "src": "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.3.tar.gz" },
    "doxygen": { "path":  thirdPartyDir+"doxygen-1.8.8", "src": "http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.8.src.tar.gz" },
    "doxytag": { "path": thirdPartyDir+"doxytag.py", "src": "https://raw.githubusercontent.com/vlfeat/vlfeat/master/docsrc/doxytag.py" },
    "cfitsio": { "path": thirdPartyDir+"cfitsio-3.36", "src": "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3360.tar.gz" } }

def Clear() :
    os.system( "rm -rf "+thirdPartyDir+"/*" )

def Main( argv ):
    usage = "Download and install third party library libraries:\n\t"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-d", u"--dry", help="Clear thirdparty dir before download and installation", action="store_true",  dest="Dry")
    parser.add_option(u"-c", u"--clear", help="Clean temporary files after download and installation", action="store_true",  dest="Clear")
    (options, args) = parser.parse_args()

    if options.Clear == True:
        Clear()
        return

    if options.Dry == True:
        Clear()

    # CPPUNIT
    libPath = libDict["cppunit"]["path"];
    libSrc = libDict["cppunit"]["src"];
    if DownloadHTTPFile( libSrc, libPath+".tar.gz" ) :
        if ExtractTarGZ( libPath+".tar.gz", libPath ) :
            os.system( "cd "+libPath+"/ ; ./configure --prefix="+thirdPartyDir+" ; make ; make install ; cd ../" )

    # GSL
    libPath = libDict["gsl"]["path"];
    libSrc = libDict["gsl"]["src"];
    if DownloadHTTPFile( libSrc, libPath+".tar.gz" ) :
        if ExtractTarGZ( libPath+".tar.gz", libPath ) :
            os.system( "cd "+libPath+"/ ; ./configure --prefix="+thirdPartyDir+" --disable-static ; make ; make install ; cd ../" )

    # FFTW
    libPath = libDict["fftw"]["path"];
    libSrc = libDict["fftw"]["src"];
    if DownloadHTTPFile( libSrc, libPath+".tar.gz" ) :
        if ExtractTarGZ( libPath+".tar.gz", libPath ) :
            os.system( "cd "+libPath+"/ ; ./configure --prefix="+thirdPartyDir+" --enable-shared; make ; make install ; cd ../" )

    # DOXYGEN
    libPath = libDict["doxygen"]["path"];
    libSrc = libDict["doxygen"]["src"];
    if DownloadHTTPFile( libSrc, libPath+".tar.gz" ) :
        if ExtractTarGZ( libPath+".tar.gz", libPath ) :
            os.system( "cd " + libPath + "/ ; ./configure --prefix " + thirdPartyDir + " ; make ; make install ; cd ../" )

    # BOOST
    libPath = libDict["boost"]["path"];
    libSrc = libDict["boost"]["src"];
    if DownloadHTTPFile( libSrc, libPath+".tar.gz" ) :
        if ExtractTarGZ( libPath+".tar.gz", libPath ) :
            os.system( "cd "+libPath+"/ ; ./bootstrap.sh --with-libraries=test,filesystem,program_options,thread,regex,python,timer,chrono --prefix="+thirdPartyDir+" ; ./b2 link=shared install ; cd ../" )

    # CFITSIO
    libPath = libDict["cfitsio"]["path"];
    libSrc = libDict["cfitsio"]["src"];
    if DownloadHTTPFile( libSrc, libPath+".tar.gz" ) :
        if ExtractTarGZ( libPath+".tar.gz", libPath ) :
            os.system( "cd "+libPath+"/ ; ./configure --enable-reentrant --prefix="+thirdPartyDir+" --enable-sse2 --enable-ssse3; make all-nofitsio; make install; cd ../" )



    # DOXYTAG
    libPath = libDict["doxytag"]["path"];
    libSrc = libDict["doxytag"]["src"];
    if DownloadHTTPFile( libSrc, libPath ) :
        os.system( "python2.7 -m compileall "+libPath )
        os.system( "mv "+thirdPartyDir+"/doxytag.pyc "+libDict["doxygen"]["path"]+"/bin/doxytag" )


Main( sys.argv )






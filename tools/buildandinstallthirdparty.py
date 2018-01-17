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
        print("File or folder: "+os.path.normpath( destPath )+" already exist, extraction skipped...")
        return False

    else :

        tfile = tarfile.open( tarPath, 'r:gz')

        print("Extracting: \n\tFrom: " + tarPath + " to: " + destPath)

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
        print("File: "+localFilePath+" already exist, download skipped...")
        return False

    # Do HTTP request
    try :
        urlfile = urllib2.urlopen( fileUrl )
    except urllib2.URLError as e:
        print("Download from: " + fileUrl + "failed.\nReason are:" + str( e.reason ))
        return

    if urlfile.info().has_key('Content-Length'):
        totalFileSize = int (urlfile.info()['Content-Length'] )
    else :
        totalFileSize = 100

    localFile = open(localFilePath, 'wb')

    print("Downloading: \n\tFrom: " + fileUrl + "\n\tTo: " + localFilePath)

    #Download file
    data_list = []
    chunk = 4096*10
    size = 0
    while 1:
        data = urlfile.read(chunk)
        if not data:
            print(" OK !")
            break
        size += len( data )
        localFile.write(data)

        progress = float( size ) / totalFileSize * 100.0
        sys.stdout.write("\r[%f%%]" %progress )

        sys.stdout.flush()

    localFile.close()

    return True

libDict = { "cppunit": { "path": "cppunit-1.12.1",
                         "src": "http://downloads.sourceforge.net/cppunit/cppunit-1.12.1.tar.gz",
                         "static": "libcppunit.a", "shared": "libcppunit.so" },
            #"boost": { "path":  "boost-1.66.0",
            #           "src": "http://downloads.sourceforge.net/project/boost/boost/1.66.0/"
            #                  "boost_1_66_0.tar.gz",
            #           "static": "libboost_chrono.a", "shared": "libboost_chrono.so.1.66.0" },
            "boost": { "path":  "boost-1.57.0",
                       "src": "http://downloads.sourceforge.net/project/boost/boost/1.57.0/"
                              "boost_1_57_0.tar.gz",
                       "static": "libboost_chrono.a", "shared": "libboost_chrono.so.1.57.0" },
            "gsl": { "path":  "gsl-2.1",
                     "src": "http://ftp.igh.cnrs.fr/pub/gnu/gsl/gsl-2.1.tar.gz",
                     "static": "libgsl.a", "shared": "libgsl.so.19.0.0" },
            "fftw": { "path":  "fftw-3.3.3",
                      "src": "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.3.tar.gz",
                      "static": "libfftw3.a", "shared": "libfftw3.so.3.3.2" },
            "doxygen": { "path":  "doxygen-1.8.8",
                         "src": "http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.8.src.tar.gz",
                         "static": "", "shared": "" },
            "doxytag": { "path": "doxytag.py",
                         "src": "https://raw.githubusercontent.com/"
                                "vlfeat/vlfeat/master/docsrc/doxytag.py",
                         "static": "", "shared": "" },
            "cfitsio": { "path": "cfitsio-3.36",
                         "src": "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/"
                                "cfitsio3360.tar.gz",
                         "static": "libcfitsio.a", "shared": "libcfitsio.so.2.3.36" }
}

def Clear() :
    os.system( "rm -rf "+thirdPartyDir+"/*" )

def _standard_build(path, third_party_dir, options):
    enable_shared = '--enable-shared --disable-static' if options.Shared else '--disable-shared --enable-static'
    os.system( "cd {path} ; ./configure --prefix={third_party_dir} {shared}; make -j{parallel}; make install ; cd ../".format(
        path=path, third_party_dir=third_party_dir, parallel=options.Parallel, shared=enable_shared) )

def _check_lib(name, shared, third_party_dir):
    filename = libDict[name]['shared' if shared else 'static']
    return os.path.exists(os.path.join(third_party_dir, 'lib', filename))


def Main( argv ):
    usage = "Download and install third party library libraries:\n\t"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-d", u"--dry", help="Clear thirdparty dir before download and installation", action="store_true",  dest="Dry")
    parser.add_option(u"-c", u"--clear", help="Clean temporary files after download and installation", action="store_true",  dest="Clear")
    parser.add_option(u"-j", u"--parallel", help="Parallel make flag", action="store", type="string", dest="Parallel", default="1")
    parser.add_option(u"-s", u"--shared", help="Build shared libs. Build static libs otherwise.", action="store_true", dest="Shared")
    (options, args) = parser.parse_args()

    if options.Clear == True:
        Clear()
        return

    if options.Dry == True:
        Clear()

    # CPPUNIT
    libPath = os.path.join(thirdPartyDir, libDict["cppunit"]["path"]);
    libSrc = libDict["cppunit"]["src"];
    if not _check_lib('cppunit', options.Shared, thirdPartyDir):
        DownloadHTTPFile( libSrc, libPath+".tar.gz" )
        ExtractTarGZ( libPath+".tar.gz", libPath )
        _standard_build(libPath, thirdPartyDir, options)

    # GSL
    libPath = os.path.join(thirdPartyDir, libDict["gsl"]["path"]);
    libSrc = libDict["gsl"]["src"];
    if not _check_lib('gsl', options.Shared, thirdPartyDir):
        DownloadHTTPFile( libSrc, libPath+".tar.gz" )
        ExtractTarGZ( libPath+".tar.gz", libPath )
        _standard_build(libPath, thirdPartyDir, options)

    # FFTW
    libPath = os.path.join(thirdPartyDir, libDict["fftw"]["path"]);
    libSrc = libDict["fftw"]["src"];
    if not _check_lib('fftw', options.Shared, thirdPartyDir):
        DownloadHTTPFile( libSrc, libPath+".tar.gz" )
        ExtractTarGZ( libPath+".tar.gz", libPath )
        _standard_build(libPath, thirdPartyDir, options)

    # DOXYGEN
    # libPath = os.path.join(thirdPartyDir, libDict["doxygen"]["path"]);
    # libSrc = libDict["doxygen"]["src"];
    # if not _check_lib('doxygen', options.Shared, thirdPartyDir):
    #     DownloadHTTPFile( libSrc, libPath+".tar.gz" )
    #     ExtractTarGZ( libPath+".tar.gz", libPath )
    #     _standard_build(libPath, thirdPartyDir, options)

    # BOOST
    libPath = os.path.join(thirdPartyDir, libDict["boost"]["path"]);
    libSrc = libDict["boost"]["src"];
    if not _check_lib('boost', options.Shared, thirdPartyDir):
        DownloadHTTPFile( libSrc, libPath+".tar.gz" )
        ExtractTarGZ( libPath+".tar.gz", libPath )
        link_type = 'shared' if options.Shared else 'static'
        os.system( "cd {path}; ./bootstrap.sh --with-libraries=test,filesystem,program_options,thread,"
                   "regex,python,timer,chrono --prefix={third_party_dir};"
                   "./b2 -j{parallel} link={link_type} install; cd ../".format(
                       path=libPath, third_party_dir=thirdPartyDir, parallel=options.Parallel,
                       link_type=link_type) )


    # CFITSIO
    libPath = os.path.join(thirdPartyDir, libDict["cfitsio"]["path"]);
    libSrc = libDict["cfitsio"]["src"];
    if not _check_lib('cfitsio', options.Shared, thirdPartyDir):
        DownloadHTTPFile( libSrc, libPath+".tar.gz" )
        ExtractTarGZ( libPath+".tar.gz", libPath )
        build_shared = 'shared' if options.Shared else 'all-nofitsio'
        os.system( "cd {path}; ./configure --enable-reentrant --prefix={third_party_dir} --enable-sse2"
                   "--enable-ssse3;make -j{parallel} {build_shared}; make install; cd ../".format(
                       path=libPath, third_party_dir=thirdPartyDir, parallel=options.Parallel,
                       build_shared=build_shared) )



    # DOXYTAG
    #libPath = os.path.join(thirdPartyDir, libDict["doxytag"]["path"]);
    #libSrc = libDict["doxytag"]["src"];
    #if DownloadHTTPFile( libSrc, libPath ) :
    #    os.system( "python2.7 -m compileall "+libPath )
    #    os.system( "mv "+thirdPartyDir+"/doxytag.pyc "+libDict["doxygen"]["path"]+"/bin/doxytag" )


Main( sys.argv )






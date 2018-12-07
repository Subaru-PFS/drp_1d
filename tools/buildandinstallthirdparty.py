#! /usr/bin/python

import urllib2
import sys
import argparse
import os
import tarfile
import shutil
from platform import platform

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

def _check_lib(name, third_party_dir, options):
    _os = platform(terse=True).lower()
    if _os == 'macos':
        ext = '.dylib' if options.shared else '.a'
    elif _os == 'windows':
        ext = '.dll' if options.shared else '.a'
    else:
        ext = '.so' if options.shared else '.a'
    filename = libDict[name]['check_file']
    return os.path.exists(os.path.join(third_party_dir, 'lib', filename + ext))

def _standard_build(path, third_party_dir, options, extra_flags=''):
    os.system( "cd {path} ; ./configure --prefix={third_party_dir} {shared}; "
               "make -j{parallel} {extra_flags} clean all; make install; cd ../".format(
        path=path, third_party_dir=third_party_dir, parallel=options.parallel,
        shared='--enable-shared' if options.shared else '--enable-static',
        extra_flags=extra_flags))

def _boost_build(path, third_party_dir, options, extra_flags=''):
    os.system( "cd {path}; ./bootstrap.sh --with-libraries=test,filesystem,program_options,thread,"
               "regex,python,timer,chrono --prefix={third_party_dir};"
               "./b2 -j{parallel} link={shared} install; cd ../".format(
                   path=path, third_party_dir=third_party_dir, parallel=options.parallel,
                   shared='shared' if options.shared else 'static') )

def _cfitsio_build(path, third_party_dir, options, extra_flags=''):
    os.system( "cd {path}; ./configure --enable-reentrant --prefix={third_party_dir} --enable-sse2 "
               "--enable-ssse3;make -j{parallel} clean {shared}; make install; cd ../".format(
                   path=path, third_party_dir=third_party_dir, parallel=options.parallel,
                   shared='shared' if options.shared else 'all-nofitsio') )

libDict = {
    #"boost": { "path":  "boost-1.66.0",
    #           "src": "http://downloads.sourceforge.net/project/boost/boost/1.66.0/"
    #                  "boost_1_66_0.tar.gz",
    #           "check_file": "libboost_chrono.a" },
    "boost": {
        "path":  "boost-1.57.0",
        "src": "http://downloads.sourceforge.net/project/boost/boost/1.57.0/"
        "boost_1_57_0.tar.gz",
        "check_file": "libboost_chrono",
        "build": _boost_build,
        "extra_flags": ''
    },
    "gsl": {
        "path":  "gsl-2.1",
        "src": "http://ftp.igh.cnrs.fr/pub/gnu/gsl/gsl-2.1.tar.gz",
        "check_file": "libgsl",
        "build": _standard_build,
        "extra_flags": ''
    },
    "fftw": {
        "path":  "fftw-3.3.3",
        "src": "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.3.tar.gz",
        "check_file": "libfftw3",
        "build": _standard_build,
        "extra_flags": "--enable-sse2 --enable-avx --enable-openmp"
    },
    "cfitsio": {
        "path": "cfitsio-3.36",
        "src": "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/"
        "cfitsio3360.tar.gz",
        "check_file": "libcfitsio",
        "build": _cfitsio_build,
        "extra_flags": ''
    }
}

def Clear() :
    os.system( "rm -rf "+thirdPartyDir+"/*" )


def Main( argv ):
    usage = "Download and install third party library libraries:\n\t"
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-d", "--dry", help="Clear thirdparty dir before download and installation",
                        action="store_true")
    parser.add_argument("-c", "--clear", help="Clean temporary files after download and installation",
                        action="store_true")
    parser.add_argument("-j", "--parallel", help="Parallel make flag", type=int,
                        default="1")
    parser.add_argument("--static", help="Build static libraries", action="store_false",
                        dest="shared", default=False)
    parser.add_argument("--shared", help="Build shared libraries", action="store_true",
                        dest="shared", default=False)
    parser.add_argument('modules', metavar='NAME', choices = libDict.keys(), nargs='*',
                        help="Modules to build.")
    args = parser.parse_args()

    if args.clear == True:
        Clear()

    if args.dry == True:
        Clear()

    for module in args.modules:
        libPath = os.path.join(thirdPartyDir, libDict[module]["path"]);
        libSrc = libDict[module]["src"];
        build_method = libDict[module]['build']
        extra_flags = libDict[module]['extra_flags']
        if not _check_lib(module, thirdPartyDir, args):
            DownloadHTTPFile(libSrc, libPath + ".tar.gz")
            ExtractTarGZ(libPath + ".tar.gz", libPath)
            build_method(libPath, thirdPartyDir, args, extra_flags)

Main( sys.argv )






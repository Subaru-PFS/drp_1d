#!/usr/bin/env python

from urllib import request, error
import sys
import argparse
import os
import tarfile
from platform import platform


def ExtractTarGZ(tarPath, destPath):
    if os.path.exists(destPath):
        print("File or folder: " + os.path.normpath(destPath) +
              " already exist, extraction skipped...")
        return False
    else:
        tfile = tarfile.open(tarPath, 'r:gz')

        print("Extracting: \n\tFrom: " + tarPath + " to: " + destPath)

        extractDir = os.path.dirname(destPath) + "/"
        tfile.extractall(extractDir)
        extractedPath = os.path.normpath(os.path.join(extractDir,
                                                      tfile.getmembers()[0].name))

        os.rename(extractedPath, destPath)

        return True


def DownloadHTTPFile(fileUrl, localFilePath):
    if localFilePath is None:
        localFilePath = os.path.basename(fileUrl)

    # Check if file already exist
    if os.path.exists(localFilePath):
        print("File: " + localFilePath + " already exist, download skipped...")
        return False

    # Do HTTP request
    try:
        urlfile = request.urlopen(fileUrl)
    except error.URLError as e:
        print("Download from: " + fileUrl + "failed.\nReason are:" + str(e.reason))
        return

    localFile = open(localFilePath, 'wb')

    print("Downloading: \n\tFrom: " + fileUrl + "\n\tTo: " + localFilePath)

    # Download file
    chunk = 4096*10
    size = 0
    while 1:
        data = urlfile.read(chunk)
        if not data:
            print(" OK !")
            break
        size += len(data)
        localFile.write(data)

    localFile.close()

    return True


def _check_lib(name, prefix, options):
    _os = platform(terse=True).lower()
    if _os == 'macos':
        ext = '.dylib' if options.shared else '.a'
    elif _os == 'windows':
        ext = '.dll' if options.shared else '.a'
    else:
        ext = '.so' if options.shared else '.a'
    filename = libDict[name]['check_file']
    return os.path.exists(os.path.join(prefix, 'lib', filename + ext))


def _standard_build(path, prefix, options, extra_flags=''):
    os.system("cd {path} ; ./configure --prefix={prefix} {shared} {extra_flags};"
              "make -j{parallel} clean all; make install".format(
                  path=path, prefix=prefix,
                  parallel=options.parallel,
                  shared='--enable-shared' if options.shared else '--enable-static',
                  extra_flags=extra_flags))


def _boost_build(path, prefix, options, extra_flags=''):
    os.system("cd {path}; ./bootstrap.sh "
              "--with-libraries=test,filesystem,program_options,thread,"
              "regex,python,timer,chrono --prefix={prefix};"
              "./b2 -j{parallel} link={shared} install; cd ../".format(
                  path=path, prefix=prefix,
                  parallel=options.parallel,
                  shared='shared' if options.shared else 'static'))


def _cfitsio_build(path, prefix, options, extra_flags=''):
    os.system("cd {path}; ./configure --enable-reentrant --prefix={prefix} "
              "--enable-sse2 --enable-ssse3 ;"
              "make -j{parallel} clean {shared}; make install; cd ../".format(
                  path=path, prefix=prefix,
                  parallel=options.parallel,
                  shared='shared' if options.shared else 'all-nofitsio'))


libDict = {
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


def Main(argv):
    if sys.version_info[0] < 3:
        print("Python version 3.X needed")
        return

    if sys.version_info[1] < 4:
        # __file__ was relative until python 3.4
        build_dir = os.path.abspath(os.path.join(os.getcwd(),
                                                 os.path.dirname(__file__),
                                                 '..', 'thirdparty'))
    else:
        build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                 '..', 'thirdparty'))

    usage = "Download and install third party library libraries:\n\t"
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("--prefix", metavar='PREFIX',
                        help="install files in PREFIX")
    parser.add_argument("-j", "--parallel",
                        help="Parallel make flag", type=int,
                        default="1")
    parser.add_argument("--static",
                        help="Build static libraries", action="store_false",
                        dest="shared", default=False)
    parser.add_argument("--shared",
                        help="Build shared libraries", action="store_true",
                        dest="shared", default=False)
    parser.add_argument('modules', metavar='NAME', choices=libDict.keys(),
                        nargs='*', help="Modules to build.")
    args = parser.parse_args()

    for module in args.modules:
        libPath = os.path.join(build_dir, libDict[module]["path"])
        libSrc = libDict[module]["src"]
        build_method = libDict[module]['build']
        extra_flags = libDict[module]['extra_flags']
        if not _check_lib(module, args.prefix, args):
            DownloadHTTPFile(libSrc, libPath + ".tar.gz")
            ExtractTarGZ(libPath + ".tar.gz", libPath)
            build_method(libPath, args.prefix, args, extra_flags)


Main(sys.argv)






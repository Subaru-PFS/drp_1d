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
        return


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

libList = []

cfitsio = "cfitsio-3.36"

libList.append(cfitsio)


def Cleanup() :
    
    for lib in libList :
        
        if os.path.exists( thirdPartyDir+lib ) :
            shutil.rmtree( thirdPartyDir+lib )
            
        if os.path.exists( thirdPartyDir+lib+".tar.gz" ) :
            os.remove( thirdPartyDir+lib+".tar.gz" )
            


def Main( argv ):

    usage = "Download and install thirdparty libraries:\n\t"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-d", u"--dry", help="Clean everything before download and installation", action="store_true",  dest="Dry")
    parser.add_option(u"-c", u"--clean", help="Clean everything after download and installation", action="store_true",  dest="Clean")
    parser.add_option(u"-C", u"--clean-only", help="Clean everything", action="store_true",  dest="CleanOnly")
    (options, args) = parser.parse_args()
    
    if options.CleanOnly == True:
        Cleanup()
        return

    if options.Dry == True:
        Cleanup()
        
    # CFITSIO
    DownloadHTTPFile( "http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3360.tar.gz", thirdPartyDir+cfitsio+".tar.gz" )
    if ExtractTarGZ( thirdPartyDir+cfitsio+".tar.gz", thirdPartyDir+cfitsio ) :
        os.system( "cd "+thirdPartyDir+cfitsio+"/ ; ./configure --enable-reentrant --prefix="+thirdPartyDir+" ; make ; make install ; cd ../" )

    if options.Clean == True:
        Cleanup()


Main( sys.argv )




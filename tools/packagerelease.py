import urllib2
import sys
import optparse
import os
import tarfile
import shutil 
import subprocess

cli_bin_name = "ez_rework_cli"
binDir = "../bin/"
buildDir = "../build/"

def CreateTargGZ( inputFiles, version, system_name, system_processor ):

    tar = tarfile.open(binDir+cli_bin_name+"_"+system_name+"_"+system_processor+"-"+version+".tar.gz", "w:gz")
    for name in inputFiles:
        if os.path.exists( name ) == False:
            print "Error: the file "+name+" does not exist, packaging aborted"
            sys.exit(1)
        tar.add(name)
    
    tar.close()



def Main( argv ):

    usage = "Build everything and package binary in a targz for distribution"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-v", u"--version", help="Specifie which version must be checkouted and builded before packaging (ex: 0.1.1) If this options is not provided, the last version will be used", action="store",  dest="Version")
    parser.add_option(u"-l", u"--list", help="List all the available version that can be packaged", action="store_true",  dest="List")
    (options, args) = parser.parse_args() 
    
    # Get list of taged versions
    tagsStr = subprocess.check_output( "git tag -l | sort -t. -k 1,1n -k 2,2n -k 3,3n -k 4,4n", shell=True )
    
    tagList = tagsStr.splitlines( )
    
    if len( tagList ) < 1:
        print "Error: No valid tags found to create a valid package"
        sys.exit(1)
        
    # List taged versions 
    if options.List == True :
        print "Available versions:"
        for tag in tagList:
            print tag + "\t"
            
        sys.exit(1)
        
    # Retrieve requested version
    version = ""
    # If no version specified, retrieve last version
    if options.Version == None:
        version = tagList[-1]
    else:
        if options.Version in tagList:
            version = options.Version
        else :
            print "Error: The given version does not match with an existing tag."
            sys.exit( 1 )
            
    print "Creating package for version: "+ version    
    
    
    system_name_file_path = "./SYSTEM_NAME"
    system_processor_file_path = "./SYSTEM_PROCESSOR"
        
    # Move to Build directory to perform every build actions
    os.chdir( buildDir )
    
    # Remove SYSTEM_* files before anything
    if os.path.exists( system_name_file_path ) :
        os.remove( system_name_file_path)
    if os.path.exists( system_processor_file_path ) :
        os.remove( system_processor_file_path )
    
    # Checkout requested version
    os.system( "git checkout tags/"+version )
    
    #Run CMake
    ret = subprocess.call( ["cmake", "."] )
    if ret != 0 :
        print "Failed to run cmake, package creation aborted"
        sys.exit( 1 )
        
    # Run Make
    ret = subprocess.call( ["make"] )
    if ret != 0 :
        print "Failed to run make, package creation aborted"
        sys.exit( 1 )
        
    # Read SYSTEM_* information created by cmake
    if os.path.exists( system_name_file_path ) :
        system_name_file = open ( system_name_file_path, "r")
        system_name = system_name_file.read().lower()
    else:
        print "File : "+system_name_file_path+" does not exist, package creation aborted"
        sys.exit( 1 )
    
    if os.path.exists( system_processor_file_path ) :    
        system_processor_file = open ( system_processor_file_path, "r")
        system_processor = system_processor_file.read().lower()
    else:
        print "File : "+system_processor_file_path+" does not exist, package creation aborted"
        sys.exit( 1 )
        
    # Move to Bin directory to perform every packaging operations
    os.chdir( "../bin/" )
    CreateTargGZ( [cli_bin_name+"_"+system_name+"_"+system_processor, "config.txt"], version, system_name, system_processor )
    
    # Checkout develop branch
    os.system( "git checkout develop" )

Main( sys.argv )

 

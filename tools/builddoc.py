import sys
import os
import optparse

def Process():
    os.chdir( "../docs/" )
    os.system( "../thirdparty/bin/doxygen ./doxygen.conf" )

def StartFromCommandLine( argv ) :	
    usage = """usage: %prog
    Run doxygen and build documentation"""
    parser = optparse.OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    
    if( len( args ) == 0 ) :
        Process( );
    else :
        print("Error: invalid argument count")
        exit()


def Main( argv ) :	
    try:
        if len( argv ) == 1 :
            StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()
    
Main( sys.argv )

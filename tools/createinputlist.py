import sys
import os
import optparse

def Process( spectrumPrefix, noisePrefix, Count, Output ) :
    f = open( Output, "w" )
	
    for i in range( 0, Count ) :
	s = spectrumPrefix+str(i)
        
        if( len( noisePrefix ) == 0 ) :
            f.write( s+".fits\n" )
        else :
	    n = noisePrefix+str(i)
            f.write( s+".fits\t"+n+".fits\n" )


def StartFromCommandLine( argv ) :	
    usage = """usage: %prog [options] spectrumPrefix noisePrefix
    ex: python ./createinputlist.py --count=1000 EZ_fits-W-TF_ EZ_fits-W-ErrStat_ """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option(u"-c", u"--count", type="int", help="Number of pair spectrum/noise",  dest="Count", default=1000 )
    parser.add_option(u"-o", u"--output", help="Output file name",  dest="Output", default="output.list")
    (options, args) = parser.parse_args()
    
    if( len( args ) == 2 ) :
        Process( args[0], args[1], options.Count, options.Output );
    elif( len( args ) == 1 ) :
        Process( args[0], "", options.Count, options.Output );
    else :
        print("Error: invalid argument count")
        exit()
    

def Main( argv ) :	
    try:
        if len( argv ) == 1 :
            #StartInteractive()
	    pass
        else:
            StartFromCommandLine( argv )
    except (KeyboardInterrupt):
        exit()
    
Main( sys.argv )

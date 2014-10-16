import sys
import pyfits
import numpy as np
import os
import optparse

magKeyName = "MAG"

def OpenFits(file2open):
	hdulist=pyfits.open(file2open)
	hdulist.info()
	return hdulist

def CloseFits(hdulist):
	hdulist.close()

def CheckInputFile( hdulist ) :
	data = hdulist[1].data
	try :
	    for i in range (len(data.field(0))):
	        z = data[i].field('Z')
	except KeyError:
	    print("Error: Input fits does not contain Z field")
      	    exit(0)

        global magKeyName

        magKeyName= "MAG"

	try :
	    for i in range (len(data.field(0))):
	        z = data[i].field( magKeyName )

	except KeyError:

            print("No MAG key present, trying with AMP key...")
            magKeyName = "AMP"

            try :
                for i in range (len(data.field(0))):
                    z = data[i].field( magKeyName )

            except KeyError:
                print("Error: Input fits does not contain MAG or AMP field")
                exit(0)



# Un fichier ascii avec 2 colonnes, 1 ligne par objet avec: Identifiant, Z et MAG
def WriteAsciiFile(hdulist,fileout):
	data = hdulist[1].data
	cols = hdulist[1].columns
	f = open(fileout, 'w')
	filename = '  id   ',cols[0].name,'   ',cols[1].name
	
	f.write('  id   '+cols[0].name+'      '+cols[1].name+'\n')
	for i in range (len(data.field(0))):
		line2print='{0}     {1}    {2}\n'.format(str(i),str(data[i].field('Z')),str(data[i].field( magKeyName )))
		f.write(line2print)
	f.close()
	return filename

# produire des spectres 1D en fits individuels correspondant a:
def WriteFitsFile(hdulist,fileout,Xcol,Ycol):
    #hdulist,fileout,Xcol,Ycol=hdu,fitsWF,'WAVE','FLUX'
	# WAVE - FLUX       #- Flux vs. wavelength  (spectre 1D simule)
	# WAVE - ERRFLUX    #- Err vs. wavelength (bruit 1D)
	# WAVE - FLUXTRUE	#- Flux_true vs. wavelength (spectre 1D de reference)
	tbdata = hdulist[1].data
	cols = hdulist[1].columns
	n = np.arange(100)
	#hdu = pyfits.PrimaryHDU(n)
	for i in range (len(cols)):
		if cols[i].name == Xcol :
			cx=cols[i]
		if cols[i].name == Ycol :
			cy=cols[i]
	
        print "Extracting Fits:"

	#for i in range (len(tbdata.field(0))):
	for i in range (len(tbdata)):
		#define the name of the output file
		destfileout=fileout+'_'+str(i)+".fits"

		print '\r'+destfileout,
                sys.stdout.flush()

		datacx=np.array(tbdata.field(Xcol)[i])
		datacy=np.array(tbdata.field(Ycol)[i])
		#print datacx
		#print datacy
		
		# define the columns to add
		col1 = pyfits.Column(name=cx.name, format='E', array=datacx*10)
		col2 = pyfits.Column(name=cy.name, format='E', array=datacy)

		# creates the coldef
		destcols=pyfits.ColDefs([col1,col2])

		# and the binary table
                #		tbhdu=pyfits.new_table(destcols)
                tbhdu=pyfits.BinTableHDU.from_columns( destcols )

		# Write in the file
		tbhdu.writeto(destfileout)

        print ""

	return destfileout



def Process( inputFile, outputDir ) :

    if os.path.exists( inputFile ) == False :
        print("Error: Input file does not exist")
	exit()

	
    hdu=OpenFits(inputFile)

    if CheckInputFile( hdu ) == False :
        print("Error: Invalid input fits format")
        exit(0)

    if len( outputDir ) == 0 :
	outputDir = os.path.splitext(inputFile)[0]+"/"
    
    if os.path.exists( outputDir ) == False :
        os.mkdir( outputDir )
    else :
        print("Error: Output dir: "+outputDir+" already exist")
	exit()

    asciifile2write=outputDir+'Z-Mag.list'
    fitsWF=outputDir+'EZ_fits-W-F'
    fitsWErrF=outputDir+'EZ_fits-W-ErrF'
    fitsWTF=outputDir+'EZ_fits-W-TF'
    fitsWErrSN=outputDir+'EZ_fits-W-SN'
    fitsWErrStat=outputDir+'EZ_fits-W-ErrStat'
    fitsWErrSys=outputDir+'EZ_fits-W-ErrSys'

    if (os.path.exists(fitsWF)):
        print("Error: Output Fits already exist")
        exit(0)

    fileascii=WriteAsciiFile(hdu,asciifile2write)
    fileWF=WriteFitsFile(hdu,fitsWF,'WAVE','FLUX')
    YES=WriteFitsFile(hdu,fitsWErrF,'WAVE','ERR')
    YES=WriteFitsFile(hdu,fitsWTF,'WAVE','FLUX_TRUE')
    YES=WriteFitsFile(hdu,fitsWErrStat,'WAVE','ERR_STAT')
    YES=WriteFitsFile(hdu,fitsWErrSys,'WAVE','ERR_SYS')

    CloseFits(hdu)

def StartFromCommandLine( argv ) :	
    usage = "usage: %prog inputFile outputDir"
    parser = optparse.OptionParser(usage=usage)
    (options, args) = parser.parse_args()
    
    if( len( args ) == 2 ) :
        Process( args[0], args[1] );
    elif( len( args ) == 1 ) :
        Process( args[0], "" );
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

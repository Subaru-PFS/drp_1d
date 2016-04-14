#include <epic/redshift/continuum/waveletsdf.h>

#include <epic/redshift/spectrum/io/fitswriter.h>
#include <epic/redshift/spectrum/io/fitsreader.h>

#include <epic/redshift/spectrum/spectrum.h>

#include <math.h>
#include <iostream>
#include <iomanip>

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;

using namespace NSEpic;
using namespace std;


CContinuumDF::CContinuumDF( std::string binPath )
{
    m_dfBinPath = binPath;
}


CContinuumDF::~CContinuumDF()
{
}

void CContinuumDF::mirror_( const CSpectrum& s, UInt32 nb, float* tab)
{
	const Float64* Ys = s.GetFluxAxis().GetSamples();
	UInt32 nn = s.GetFluxAxis().GetSamplesCount();
	UInt32 nall = nn+nb*2;
	UInt32 k;
	for (Int32 j=0; j<nb; j++) {
		k = nb-j;
		tab [j] = Ys[k];
	}

	for (Int32 j=nb; j<nall-nb; j++) {
		tab [j] = Ys[j-nb];
	}

	for (Int32 j=nall-nb; j<nall; j++) {
		k = (nall-j)+nn-nb-1;
		tab [j] = Ys[k];
	}

}


std::string CContinuumDF::temporaryPath ( const CSpectrum& s, UInt32 nall)
{
	UInt32 nb   = round(s.GetFluxAxis().GetSamplesCount()/10);

	extendedData_ = new float[nall];
	mirror_( s, nb,extendedData_ );

	bfs::path temporaPath0 = bfs::unique_path();
	temporaPath0="./tempoData"/temporaPath0;
	bfs::create_directories(temporaPath0);
	bfs::path tm0 = temporaPath0/"extendedData.fits";
	const char* filePathTEST = tm0.c_str();


	fitsfile *fptr = NULL;
	Int32 status = 0;
	Int32 hdunum=0;
	Bool retv = true;
	long naxis = 1;
	long naxes[1] = {nall};
	if( bfs::exists( filePathTEST ) ) {	bfs::remove( filePathTEST );}
	Int32 nocrea = fits_create_file(&fptr, filePathTEST, &status);
	Int32 noimg  = fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
	Int32 nowrite = fits_write_img(fptr, TFLOAT, 1, nall, extendedData_, &status);
	Int32 noclose = fits_close_file(fptr, &status);

	return filePathTEST;
}



Bool CContinuumDF::RemoveContinuum ( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis)
{
	Int32 nb   = round(s.GetFluxAxis().GetSamplesCount()/10);
	Int32 nn   = s.GetFluxAxis().GetSamplesCount();
	Int32 nall = nn+nb*2;

	Int32 nscales = s.GetDecompScales();
	std::string decomp_pmt    = boost::lexical_cast<std::string>(nscales);
	std::string decomp_spline = boost::lexical_cast<std::string>(nscales+1);


	/* *******************************************************
	*            STEP 0 :    SET TEMPORARY PATHS
	* ****************************************************** */
	bfs::path temporaPath = bfs::unique_path();
	temporaPath="./tempoData"/temporaPath;
	bfs::create_directories(temporaPath);

	bfs::path tm = temporaPath/"resPMT.fits";
	//if( bfs::exists( tm ) ) {	bfs::remove( tm );}
	std::string outputFilePMT = tm.c_str();

	bfs::path tm2 = temporaPath/"resSpline.fits";
	//if( bfs::exists( tm2 ) ) {	bfs::remove( tm2 );}
	std::string outputFileSpline = tm2.c_str();
	//std::cout<< "bfsPath() = "<< bfs::system_complete(outputFileSpline)<<std::endl;


	/* *******************************************************
	*            STEP 1 :    APPLY PMT TRANSFORM
	* ****************************************************** */
	std::string inputFile = temporaryPath(s, nall);
	//std::string filePath = "./extern/mr1d_filter_modified";
    bfs::path binPath = m_dfBinPath;
    binPath.append("mr1d_filter_modified");
    std::string filePath = binPath.string();
	std::string params = "-n"+decomp_pmt+" -t11 -f3 -i20 -k -K -s8 -P ";
	std::string cmm= filePath+" "+params+" "+inputFile+" "+outputFilePMT;
	std::string command = "exec "+ cmm;
	system(command.c_str());


	/* *******************************************************
	*            STEP 2 :    APPLY B3SPLINE TRANSFORM
	* ****************************************************** */
	std::string inputFileSpline = outputFilePMT;
    bfs::path binPathSpline = m_dfBinPath;
    binPathSpline.append("mr1d_trans");
    std::string filePathSpline = binPathSpline.string();
	std::string paramsSpline = "-n"+decomp_spline+" -t3 ";
	std::string cmmspline= filePathSpline+" "+paramsSpline+" "+inputFileSpline+" "+outputFileSpline;
	std::string commandspline = "exec "+ cmmspline;
	system(commandspline.c_str());


	/* *******************************************************
	 *           READ PMT-B3SPLINE RESULT
	 * ****************************************************** */
	fitsfile *fptr2 = NULL;
	Int32 status = 0;
	Int32 hdunum=0;
	Int32 nullval = 0;
	Int32 anynul = 0;

	float array[nscales+1][nall];
	Int32 noopen = fits_open_file(&fptr2, outputFileSpline.c_str(), READONLY, &status);
	if ( noopen )
	{
		return false;
	}
	Int32 noread  = fits_read_img(fptr2, TFLOAT, 1, nall*(nscales+1), &nullval, array, &anynul, &status);
	if ( noread )
	{
		return false;
	}

	estimatedBaseline_extd =new float[nall];
	estimatedBaseline_extd= array[nscales];

	estimatedBaseline =new float[nn];

	CSpectrum s_baseline;
	CSpectrumFluxAxis& baseline = s_baseline.GetFluxAxis();
	CSpectrumAxis& wavAxis = s_baseline.GetSpectralAxis();
	wavAxis =  s.GetSpectralAxis();

	baseline.SetSize(nn);
	CSpectrumFluxAxis fluxAxis = s.GetFluxAxis();

	noContinuumFluxAxis.SetSize(nn);
	for(Int32 j=0;j<nn;j++)
	{
		estimatedBaseline[j]        = estimatedBaseline_extd[j+nb];
		baseline[j]                       =  estimatedBaseline[j]  ;
		noContinuumFluxAxis[j] = fluxAxis[j]-estimatedBaseline[j];
	}

	Float64* noContinuumFluxAxisError = noContinuumFluxAxis.GetError();
	const Float64* fluxAxisError = fluxAxis.GetError();
	for(Int32 j=0;j<nn;j++)
	{
			noContinuumFluxAxisError[j] = fluxAxisError[j];
	}


	/* *******************************************************
	*            DEPRECATED: Local save
	* ******************************************************
	if (0)
	{
		bfs::path pp = boost::filesystem::current_path()/"/../tempoData/baselineDF/";
		std::string name = s.GetName()+"_DF.fits";
		pp=pp/name;
		CSpectrumIOFitsWriter writer;
		pp = bfs::system_complete(pp);
		Bool retVal = writer.Write(pp.c_str(), s_baseline );
	}
	*/


	/* *******************************************************
	*            REMOVE TEMPORARY PATHS
	* ****************************************************** */
	bfs::remove_all(bfs::system_complete(temporaPath));    							   // remove pmt & b3spline temporary Resuts
	bfs::remove_all(bfs::system_complete(bfs::path(inputFile).parent_path()));  // remove extended data


    return true;

}



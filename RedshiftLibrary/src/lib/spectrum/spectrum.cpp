#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>

#include <RedshiftLibrary/continuum/waveletsdf.h>
#include <RedshiftLibrary/continuum/median.h>
#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>

#include <RedshiftLibrary/log/log.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/debug/assert.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

#include <gsl/gsl_fit.h>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

using namespace NSEpic;
using namespace std;

CSpectrum::CSpectrum():
    m_estimationMethod(""),
    m_medianWindowSize(-1),
    m_nbScales(-1),
    m_dfBinPath(""),
    m_Name("")
{

}

CSpectrum::CSpectrum(const CSpectrum& other, TFloat64List mask):
    m_estimationMethod(other.m_estimationMethod),
    m_dfBinPath(other.m_dfBinPath),
    m_medianWindowSize(other.m_medianWindowSize),
    m_nbScales(other.m_nbScales),
    m_Name(other.m_Name),
    m_spcType(other.m_spcType),
    m_SpectralAxis(UInt32(0), other.m_SpectralAxis.IsInLogScale())
{
    const CSpectrumSpectralAxis & otherSpectral = other.m_SpectralAxis;
    const CSpectrumFluxAxis & otherFlux = other.GetFluxAxis();
    const TFloat64List& error = otherFlux.GetError();

    TAxisSampleList& SpectralVector = m_SpectralAxis.GetSamplesVector();
    TAxisSampleList& FluxVector = m_FluxAxis->GetSamplesVector();
    TFloat64List& ErrorVector = m_FluxAxis->GetError();

    const UInt32 otherSpectralSize = otherSpectral.GetSamplesCount();
    const UInt32 otherFluxSize = otherFlux.GetSamplesCount();
    UInt32 minsize = min((UInt32)mask.size(), otherSpectralSize);
    minsize = min(minsize, otherFluxSize);
    for(Int32 i=0; i<minsize; i++){
        if(mask[i]!=0){
            SpectralVector.push_back(otherSpectral[i]);
            FluxVector.push_back(otherFlux[i]);
            if( !error.empty() ){
                ErrorVector.push_back(error[i]);
            }
        }
    }
}

CSpectrum::CSpectrum(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis) :
    m_SpectralAxis(spectralAxis),
    m_estimationMethod(""),
    m_medianWindowSize(-1),
    m_nbScales(-1),
    m_dfBinPath(""),
    m_Name("") 
{
    *m_FluxAxis = fluxAxis;
}

//assignment constructor
CSpectrum::CSpectrum(const CSpectrum& other):
    m_estimationMethod(other.m_estimationMethod),
    m_dfBinPath(other.m_dfBinPath),
    m_medianWindowSize(other.m_medianWindowSize),
    m_nbScales(other.m_nbScales),
    m_SpectralAxis(other.m_SpectralAxis),
    m_spcType(other.m_spcType),
    m_Name(other.m_Name)
{
    // need to switch m_fluxAxis to correct buffer
    // here it points by default to raw

}

CSpectrum::~CSpectrum()
{

}

//copy constructor
CSpectrum& CSpectrum::operator=(const CSpectrum& other)
{
    m_SpectralAxis = other.m_SpectralAxis;
    *m_FluxAxis = other.GetFluxAxis();

    m_estimationMethod = other.m_estimationMethod;
    m_dfBinPath = other.m_dfBinPath;
    m_medianWindowSize = other.m_medianWindowSize;
    m_nbScales = other.m_nbScales;
    m_Name = other.m_Name;

    return *this;
}

/**
 * below should be calculated in the case of precomputedfinegrid
 */
Bool CSpectrum::RebinFineGrid() const
{
  // Precalculate a fine grid template to be used for the 'closest value' rebin method
  Int32 n = m_FluxAxis->GetSamplesCount();
  if(!n)
    return false;

  Float64 lmin = m_SpectralAxis[0];//template wavelength never starts at 0
  Float64 lmax = m_SpectralAxis[n - 1];
  Int32 nTgt = (lmax - lmin) / m_dLambdaFineGrid + 2.0 / m_dLambdaFineGrid;

  m_pfgFlux.resize(nTgt);

  const TAxisSampleList Ysrc = m_FluxAxis->GetSamplesVector();
  const TAxisSampleList Xsrc = m_SpectralAxis.GetSamplesVector();

  //Initialize and allocate the gsl objects
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spline, Xsrc.data(), Ysrc.data(), n);
  gsl_interp_accel* accelerator = gsl_interp_accel_alloc();

  Int32 k = 0;
  Float64 x = 0.0;
  for (k = 0; k < nTgt; k++) {
    x = lmin + k * m_dLambdaFineGrid;
    if (x < m_SpectralAxis[0] || x > m_SpectralAxis[n - 1]) {
      m_pfgFlux[k] = 0.0;
    } else {
      m_pfgFlux[k] = gsl_spline_eval(spline, x, accelerator);
    }
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(accelerator);
  m_FineGridInterpolated = true;
  return true;
}

Bool CSpectrum::RemoveContinuum( CContinuum& remover ) const
{
    return remover.RemoveContinuum( *this, m_WithoutContinuumFluxAxis );
}

/**
 * Estimate the continuum component
 */
void CSpectrum::EstimateContinuum() const
{
    Log.LogInfo( "Continuum estimation on input spectrum: using %s", m_estimationMethod.c_str() );
    if( m_estimationMethod == "IrregularSamplingMedian" )
    {
        CContinuumIrregularSamplingMedian continuum;
        continuum.SetMedianKernelWidth( m_medianWindowSize );
        continuum.SetMeanKernelWidth( m_medianWindowSize );
        RemoveContinuum( continuum );
        Log.LogInfo( "Continuum estimation - medianKernelWidth = %.2f", m_medianWindowSize );
    }else if( m_estimationMethod == "Median" )
    {
        CContinuumMedian continuum;
        continuum.SetMedianKernelWidth( m_medianWindowSize );
        RemoveContinuum( continuum );
        Log.LogInfo( "Continuum estimation - medianKernelWidth = %.2f", m_medianWindowSize );
    }else if( m_estimationMethod == "waveletsDF" )
    {
        CContinuumDF continuum(m_dfBinPath);
        Bool ret = RemoveContinuum( continuum );
        if( !ret ) //doesn't seem to work. TODO: check that the df errors lead to a ret=false value
        {
          Log.LogError("Failed to apply continuum substraction for spectrum: %s", this->GetName().c_str());
          throw std::runtime_error("Failed to apply continuum substraction");
        }
    }else if( m_estimationMethod == "raw" )
    {
        Int32 nbSamples = this->GetSampleCount();
        m_WithoutContinuumFluxAxis.SetSize(nbSamples);
        for(Int32 k=0; k<nbSamples; k++)
        {
            m_WithoutContinuumFluxAxis[k] = 0.0;
        }
    }else if( m_estimationMethod == "zero" )
    {
        m_WithoutContinuumFluxAxis = m_RawFluxAxis;
    }

    m_nameBaseline = m_method2baseline[m_estimationMethod];
    Log.LogInfo("===============================================");

    // Fill m_ContinuumFluxAxis
    CSpectrumFluxAxis continuumFluxAxis = m_RawFluxAxis;
    continuumFluxAxis.Subtract( m_WithoutContinuumFluxAxis );
    m_ContinuumFluxAxis = continuumFluxAxis;
}

std::string CSpectrum::GetBaseline() const
{
    return m_nameBaseline;
}

/**
 * Invert the flux axis
 */
Bool CSpectrum::InvertFlux()
{
    m_FluxAxis->Invert();
    return true;
}

/**
 * Convert the spectral axis to a neperian logarithm scale
 */
Bool CSpectrum::ConvertToLogScale()
{
    return m_SpectralAxis.ConvertToLogScale();
}

/**
 * Convert the spectral axis to a linear scale
 */
Bool CSpectrum::ConvertToLinearScale()
{
    return m_SpectralAxis.ConvertToLinearScale();
}

Float64 CSpectrum::GetResolution() const
{
    return m_SpectralAxis.GetResolution();
}

Float64 CSpectrum::GetMeanResolution() const
{
    return m_SpectralAxis.GetMeanResolution();
}

/**
 * Return the lambda range of the entire spectrum.
 * Range is always expressed in linear scale NOT in log scale even if the underlying spcetrum is in log scale
 */
TLambdaRange CSpectrum::GetLambdaRange() const
{
    return m_SpectralAxis.GetLambdaRange();
}

bool CSpectrum::GetMeanAndStdFluxInRange(TFloat64Range wlRange, Float64& mean, Float64 &std) const
{
    //wlrange should be totally included in the spectrum lambdarange
    if(wlRange.GetBegin()<m_SpectralAxis.GetLambdaRange().GetBegin())
    {
        return false;
    }
    if(wlRange.GetEnd()>m_SpectralAxis.GetLambdaRange().GetEnd())
    {
        return false;
    }

    CMask mask;
    m_SpectralAxis.GetMask( wlRange, mask );
    const TFloat64List& error = m_FluxAxis->GetError();
    Float64 _Mean = 0.0;
    Float64 _SDev = 0.0;
    m_FluxAxis->ComputeMeanAndSDev(mask, _Mean, _SDev, error);

    mean = _Mean;
    std = _SDev;
    return true;
}

bool CSpectrum::GetLinearRegInRange(TFloat64Range wlRange, Float64 &a, Float64 &b) const
{
    //wlrange should be totally included in the spectrum lambdarange
    if(wlRange.GetBegin()<m_SpectralAxis.GetLambdaRange().GetBegin())
    {
        return false;
    }
    if(wlRange.GetEnd()>m_SpectralAxis.GetLambdaRange().GetEnd())
    {
        return false;
    }

    const TFloat64List& error = m_FluxAxis->GetError();

    TInt32Range iRange = m_SpectralAxis.GetIndexesAtWaveLengthRange(wlRange);
    Int32 n = iRange.GetEnd()-iRange.GetBegin()+1;
    Float64* x = (Float64*)malloc(n*sizeof(Float64));
    Float64* y = (Float64*)malloc(n*sizeof(Float64));
    Float64* w = (Float64*)malloc(n*sizeof(Float64));

    for(Int32 k=0; k<n; k++)
    {
        Int32 ik=k+iRange.GetBegin();
        w[k] = 1.0/(error[ik]*error[ik]);
        x[k] = m_SpectralAxis[ik];
        y[k] = m_FluxAxis->GetSamples()[ik];
    }

    double c0, c1, cov00, cov01, cov11, chisq;
    gsl_fit_wlinear (x, 1, w, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

    a = c1;
    b = c0;
    //todo: use chisq (weighted sum of squares of the residuals) for fitting quality,
    free(x);
    free(y);
    free(w);
    return true;
}

const std::string CSpectrum::GetName() const
{
    return m_Name;
}

void CSpectrum::SetName(const char* name)
{
    m_Name = name;
}

const Int32 CSpectrum::GetType() const
{
    return m_spcType;
}

void CSpectrum::SetType(const Int32 type)
{
    m_spcType = type;
}

const Bool CSpectrum::checkFlux( Float64 flux, Int32 index ) const
{
    //Log.LogDebug("    CSpectrum::checkFlux - Found flux value (=%e) at index=%d", flux, index);
    Bool validValue = true;
    if( std::isnan(flux) ){
        validValue = false;
        //Log.LogDebug("    CSpectrum::checkFlux - Found nan flux value (=%e) at index=%d", flux, index);
    }
    if( std::isinf(flux) ){
        validValue = false;
        //Log.LogDebug("    CSpectrum::checkFlux - Found inf flux value (=%e) at index=%d", flux, index);
    }
    if( flux != flux ){
        validValue = false;
        //Log.LogDebug("    CSpectrum::checkFlux - Found invalid flux value (=%e) at index=%d", flux, index);
    }
    return validValue;
}

const Bool CSpectrum::checkNoise( Float64 error, Int32 index ) const
{
    Log.LogDebug("    CSpectrum::checkNoise - Found noise value (=%e) at index=%d", error, index);
    Bool validValue = true;
    if( error < DBL_MIN ){
        //check if noise is below minimum normalized positive value of double
        validValue = false;
        Log.LogDebug("    CSpectrum::checkNoise - Found subnormal noise value (=%e) at index=%d", error, index);
    }
    if( std::isnan(error) ){
        validValue = false;
        Log.LogDebug("    CSpectrum::checkNoise - Found nan noise value (=%e) at index=%d", error, index);
    }
    if( std::isinf(error) ){
        validValue = false;
        Log.LogDebug("    CSpectrum::checkNoise - Found inf noise value (=%e) at index=%d", error, index);
    }
    if( error != error ){
        validValue = false;
        Log.LogDebug("    CSpectrum::checkNoise - Found invalid noise value (=%e) at index=%d", error, index);
    }
    return validValue;
}

const Bool CSpectrum::IsFluxValid( Float64 LambdaMin, Float64 LambdaMax ) const
{
    Bool allzero = true;
    Bool invalidValue = false;
    Int32 nInvalid = 0;

    const Float64* flux = m_FluxAxis->GetSamples();
    if (LambdaMin < m_SpectralAxis[0] || LambdaMax > m_SpectralAxis[m_SpectralAxis.GetSamplesCount()-1]){
        return false;
    }
    else
    {
        Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
        Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
        Log.LogDetail( "CSpectrum::IsFluxValid - checking on the configured lambdarange = (%f, %f)", LambdaMin, LambdaMax );
        Log.LogDetail( "CSpectrum::IsFluxValid - checking on the true observed spectral axis lambdarange = (%f, %f)", m_SpectralAxis[iMin], m_SpectralAxis[iMax] );
        for(Int32 i=iMin; i<iMax; i++){
            //check flux
            Bool validSample = checkFlux(flux[i], i);

            if(!validSample){
                invalidValue = true;
                nInvalid++;
            }

            //all zero check
            if( flux[i] != 0.0 ){
                allzero = false;
                //Log.LogDebug("    CSpectrum::IsFluxValid - Found non zero and valid flux value (=%e) at index=%d", i, flux[i]);
            }
        }
        Bool valid = !invalidValue && !allzero;
        if(nInvalid>0)
        {
            Log.LogDetail("    CSpectrum::IsFluxValid - Found %d invalid flux samples", nInvalid);
        }
        return valid;
    }
}

const Bool CSpectrum::IsNoiseValid( Float64 LambdaMin, Float64 LambdaMax ) const
{
    Bool valid = true;
    Int32 nInvalid = 0;

    const TFloat64List& error = m_FluxAxis->GetError();
    if (LambdaMin < m_SpectralAxis[0] || LambdaMax > m_SpectralAxis[m_SpectralAxis.GetSamplesCount()-1]){
        return false;
    }
    else
    {
        Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
        Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
        Log.LogDetail( "CSpectrum::IsNoiseValid - checking on the configured lambdarange = (%f, %f)", LambdaMin, LambdaMax );
        Log.LogDetail( "CSpectrum::IsNoiseValid - checking on the true observed spectral axis lambdarange = (%f, %f)", m_SpectralAxis[iMin], m_SpectralAxis[iMax] );
        //Log.LogDebug("    CSpectrum::IsNoiseValid - debug - iMin=%d and wmin=%f, iMax=%d and wmax=%f", iMin, m_SpectralAxis[iMin], iMax, m_SpectralAxis[iMax]);
        for(Int32 i=iMin; i<iMax; i++){
            //check noise
            Bool validSample = checkNoise(error[i], i);

            if(!validSample){
                valid = false;
                nInvalid++;
            }
        }
        if(nInvalid>0)
        {
            Log.LogDetail("    CSpectrum::IsNoiseValid - Found %d invalid noise samples", nInvalid);
        }
        return valid;
    }
}

Bool CSpectrum::correctSpectrum( Float64 LambdaMin, Float64 LambdaMax, Float64 coeffCorr )
{
    Bool corrected = false;
    Int32 nCorrected = 0;

    TFloat64List& error = m_FluxAxis->GetError();
    Float64 *flux = m_FluxAxis->GetSamples();

    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    //Log.LogDebug("    CSpectrum::correctSpectrum - debug - iMin=%d and wmin=%f, iMax=%d and wmax=%f", iMin, m_SpectralAxis[iMin], iMax, m_SpectralAxis[iMax]);

    Float64 maxNoise = -DBL_MAX;
    Float64 minFlux = DBL_MAX;
    for(Int32 i=iMin; i<iMax; i++){
        //Log.LogDebug("    CSpectrum::correctSpectrum - debug - RAW sample for maxFlux/minNoise. Found err(=%f) and flux=%f", error[i], flux[i]);
        //check noise & flux
	if (!checkNoise(error[i], i) || !checkFlux(flux[i], i))
	    continue;

        //Log.LogDebug("    CSpectrum::correctSpectrum - debug - valid sample for maxFlux/minNoise. Found err(=%f) and flux=%f", error[i], flux[i]);
        if(error[i] > maxNoise)
        {
            maxNoise = error[i];
        }
        if(std::abs(flux[i]) < std::abs(minFlux))
        {
            minFlux = abs(flux[i]);
        }
    }
    if(minFlux==DBL_MAX)
    {
        Log.LogError("    CSpectrum::correctSpectrum - unable to set minFlux value (=%e). Setting it to 0.", minFlux);
        minFlux = 0.0;
    }
    if(maxNoise==-DBL_MAX)
    {
        Log.LogError("    CSpectrum::correctSpectrum - unable to set maxNoise value.");
        throw std::runtime_error("    CSpectrum::correctSpectrum - impossible to correct input spectrum.");
    }

    for(Int32 i=iMin; i<iMax; i++){
        //check noise & flux
        Bool validSample = checkNoise(error[i], i) && checkFlux(flux[i], i);

        if(!validSample){
            error[i] = maxNoise*coeffCorr;
            flux[i] = minFlux/coeffCorr;
            corrected=true;
            nCorrected++;
        }
    }
    if(nCorrected>0)
    {
        Log.LogInfo("    CSpectrum::correctSpectrum - Corrected %d invalid samples with coeff (=%f), minFlux=%e, maxNoise=%e", nCorrected, coeffCorr, minFlux, maxNoise);
    }
    return corrected;
}

const std::string& CSpectrum::GetFullPath() const
{
      return m_FullPath;
}

const Int32 CSpectrum::GetDecompScales() const
{
      return m_nbScales;
}

const Float64 CSpectrum::GetMedianWinsize() const
{
      return m_medianWindowSize;
}

const std::string CSpectrum::GetContinuumEstimationMethod() const
{
      return m_estimationMethod;
}

const std::string CSpectrum::GetWaveletsDFBinPath() const
{
      return m_dfBinPath;
}

void CSpectrum::SetFullPath(const char* nameP)
{
    m_FullPath = nameP;
}

void CSpectrum::SetDecompScales( Int32 decompScales )
{
    m_nbScales = decompScales;
}

void CSpectrum::SetMedianWinsize( Float64 winsize )
{
    m_medianWindowSize = winsize;
}

void CSpectrum::SetContinuumEstimationMethod( std::string method )
{
    m_estimationMethod = method;
}

void CSpectrum::SetWaveletsDFBinPath(std::string binPath)
{
    m_dfBinPath = binPath;
}

void CSpectrum::LoadSpectrum(const char* spectrumFilePath, const char* noiseFilePath)
{
  std::string spcName="";
  std::string noiseName="";

  CSpectrumIOGenericReader reader;

  if(spectrumFilePath!=NULL)
      spcName = bfs::path(spectrumFilePath).stem().string() ;

  if(noiseFilePath!=NULL)
      noiseName = bfs::path(noiseFilePath).stem().string() ;

  SetName(spcName.c_str());
  SetFullPath(bfs::path( spectrumFilePath ).string().c_str());

  reader.Read(spectrumFilePath, *this);
  Log.LogInfo("Successfully loaded input spectrum file: (%s), samples : %d",
	      spcName.c_str(), GetSampleCount() );

  // add noise if any or add flat noise
  if( noiseFilePath == NULL )
    {
      CNoiseFlat noise;
      noise.SetStatErrorLevel( 1.0 );
      noise.AddNoise(*this);
    }
  else
    {
      CNoiseFromFile noise;
      noise.SetNoiseFilePath(noiseFilePath, reader);
      noise.AddNoise(*this);
      Log.LogInfo("Successfully loaded input noise file:    (%s)", noiseName.c_str() );
    }
}

void CSpectrum::InitPrecomputeFineGrid() const
{
    m_FineGridInterpolated = false;
}

///
/// * This rebin method targets processing speed:
/// - it uses already allocated rebinedFluxAxis, rebinedSpectralAxis and rebinedMask
/// - opt_interp = 'lin' : linear interpolation is performed by default
/// - opt_interp = 'precomputedfinegrid' : nearest grid point interpolation is performed using m_pfgFlux which is the precomputed fine grid
/// - opt_interp = 'spline' : GSL/spline interpolation is performed (TODO - not tested)
/// - opt_interp = 'ngp' : nearest grid point is performed (TODO - not tested)
/**
 * targetSpectralAxis should be expressed in same frame as source SpetralAxis
*/
Bool CSpectrum::Rebin( const TFloat64Range& range, const CSpectrumSpectralAxis& targetSpectralAxis,
                       CSpectrum& rebinedSpectrum, CMask& rebinedMask, const std::string opt_interp, const std::string opt_error_interp ) const
{
    
    if( m_SpectralAxis.GetSamplesCount() != m_FluxAxis->GetSamplesCount() )
    {
        Log.LogError("Problem samplecountsize betwees spectral axis and flux axis");
        return false;
    }
    
    if( opt_interp=="precomputedfinegrid" && m_FineGridInterpolated == false )
    {
        RebinFineGrid();
    }

    if(m_pfgFlux.size()==0 && opt_interp=="precomputedfinegrid" && m_FineGridInterpolated == true)
    {
        Log.LogError("Problem buffer couldnt be computed\n" );
        return false;
    }

    //the spectral axis should be in the same scale
    TFloat64Range logIntersectedLambdaRange( log( range.GetBegin() ), log( range.GetEnd() ) );
    TFloat64Range currentRange = logIntersectedLambdaRange;
    if(m_SpectralAxis.IsInLinearScale() != targetSpectralAxis.IsInLinearScale() ){
        Log.LogError("Problem spectral axis and target spectral axis are not in same scale\n");
        return false;

    }
    if(m_SpectralAxis.IsInLinearScale()){
        currentRange = range;
    }
    UInt32 s = targetSpectralAxis.GetSamplesCount();

    rebinedSpectrum.m_SpectralAxis = targetSpectralAxis; // copy (necessary)

    CSpectrumFluxAxis& rebinedFluxAxis = rebinedSpectrum.GetFluxAxis();
    rebinedFluxAxis.SetSize(s);  // does not re-allocate if already allocated

    rebinedMask.SetSize(s);

    const TAxisSampleList& Xsrc = m_SpectralAxis.GetSamplesVector();
    const TAxisSampleList& Ysrc = m_FluxAxis->GetSamplesVector();
    const TAxisSampleList& Xtgt = targetSpectralAxis.GetSamplesVector();
    TAxisSampleList&       Yrebin = rebinedFluxAxis.GetSamplesVector();
    const TFloat64List&    Error = m_FluxAxis->GetError();
    TFloat64List&          ErrorRebin = rebinedFluxAxis.GetError();

    // Move cursors up to lambda range start
    Int32 j = 0;
    while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] < currentRange.GetBegin() )
    {
        rebinedMask[j] = 0;
        Yrebin[j] = 0.0;
        if (opt_error_interp == "rebin")
            ErrorRebin[j] = INFINITY;
        j++;
    }

    if(opt_interp=="lin"){
        //Default linear interp.
        Int32 k = 0;
        // For each sample in the valid lambda range interval.
        while( k<m_SpectralAxis.GetSamplesCount()-1 && Xsrc[k] <= currentRange.GetEnd() )
        {           
            // For each sample in the target spectrum that are in between two continous source sample
            while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= Xsrc[k+1] )
            {
                // perform linear interpolation of the flux
                Float64 xSrcStep = ( Xsrc[k+1] - Xsrc[k] );
                Float64 t = ( Xtgt[j] - Xsrc[k] ) / xSrcStep;
                Yrebin[j] = Ysrc[k] + ( Ysrc[k+1] - Ysrc[k] ) * t;
                rebinedMask[j] = 1;

                if (opt_error_interp != "no")
                {
                    if (opt_error_interp == "rebin" )
                        ErrorRebin[j] = Error[k] + ( Error[k+1] - Error[k] ) * t;
                    else if (opt_error_interp == "rebinVariance")
                    {
                        ErrorRebin[j] = sqrt(Error[k]*Error[k]*(1-t)*(1-t) + Error[k+1]*Error[k+1] * t*t);
                        //*
                        Float64 xDestStep, xStepCompensation;
                        if(j<targetSpectralAxis.GetSamplesCount()-1)
                        {
                            xDestStep = Xtgt[j+1]-Xtgt[j];
                            xStepCompensation = xSrcStep/xDestStep;
                        }else if(j>0){
                            xDestStep = Xtgt[j]-Xtgt[j-1];
                            xStepCompensation = xSrcStep/xDestStep;
                        }else{
                            xStepCompensation = 1.0;
                        }
                        ErrorRebin[j] *= sqrt(xStepCompensation);
                    }
                }
                j++;
            }

            k++;
        }
    }else if(opt_interp=="precomputedfinegrid"){
        // Precomputed FINE GRID nearest sample, 20150801
        Int32 k = 0;
        Float64 lmin = m_SpectralAxis[0];
       // For each sample in the target spectrum
        while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
        {
            k = int( (Xtgt[j] - lmin)/m_dLambdaFineGrid + 0.5 );
            Yrebin[j] = m_pfgFlux[k];
            rebinedMask[j] = 1;

            // note: error rebin not implemented for precomputedfinegrid
            if (opt_error_interp != "no")
                return false;

            j++;

        }
    }else if(opt_interp=="spline"){
        // GSL method spline
        //Initialize and allocate the gsl objects
        Int32 n = m_SpectralAxis.GetSamplesCount();
        gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, n);
        gsl_spline_init (spline, Xsrc.data(), Ysrc.data(), n);
        gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();

        // For each sample in the valid lambda range interval.
        while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
        {
            Yrebin[j] = gsl_spline_eval (spline, Xtgt[j], accelerator);
            rebinedMask[j] = 1;

            // note: error rebin not implemented for spline interp
            if (opt_error_interp != "no")
                return false;

            j++;
        }
    }else if(opt_interp=="ngp"){
        //nearest sample, lookup
        Int32 k = 0; 
        Int32 kprev = 0;
        Int32 n = m_SpectralAxis.GetSamplesCount();
        while( j<targetSpectralAxis.GetSamplesCount() && Xtgt[j] <= currentRange.GetEnd() )
        {
            k = gsl_interp_bsearch (Xsrc.data(), Xtgt[j], kprev, n);
            kprev = k;
            // closest value
            Yrebin[j] = Ysrc[k];

            if (opt_error_interp != "no")
            {
                ErrorRebin[j] = Error[k];

                if (opt_error_interp == "rebinVariance")
                {
                     Float64 xSrcStep, xDestStep, xStepCompensation;
                     if(j<targetSpectralAxis.GetSamplesCount()-1)
                     {
                         xDestStep = Xtgt[j+1]-Xtgt[j];
                         xStepCompensation = xSrcStep/xDestStep;
                     }else if(j>0){
                         xDestStep = Xtgt[j]-Xtgt[j-1];
                         xStepCompensation = xSrcStep/xDestStep;
                     }else{
                         xStepCompensation = 1.0;
                     }
                     ErrorRebin[j] *= sqrt(xStepCompensation);
                }
            }

            rebinedMask[j] = 1;
            j++;
        }
    }

    while( j < targetSpectralAxis.GetSamplesCount() )
    {
        rebinedMask[j] = 0;
        Yrebin[j] = 0.0;
        if (opt_error_interp == "rebin" )
            ErrorRebin[j] = INFINITY;
        j++;
    }

    return true;
}

void CSpectrum::ScaleFluxAxis(Float64 scale){
    *m_FluxAxis *= scale;
}

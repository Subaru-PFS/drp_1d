#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>

#include <RedshiftLibrary/log/log.h>

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

CSpectrum::CSpectrum()
{
    m_estimationMethod = "";
    m_medianWindowSize = -1;
    m_nbScales = -1;
    m_dfBinPath = "";

}


CSpectrum::CSpectrum(const CSpectrum& other, TFloat64List mask)
{
    TFloat64List tmpFlux;
    TFloat64List tmpError;
    TFloat64List tmpWave;

    CSpectrumSpectralAxis otherSpectral = other.GetSpectralAxis();
    CSpectrumFluxAxis otherFlux = other.GetFluxAxis();

    const TFloat64List& error = otherFlux.GetError();

    for(Int32 i=0; i<mask.size(); i++){
        if(mask[i]!=0){
            tmpWave.push_back(otherSpectral[i]);
            tmpFlux.push_back(otherFlux[i]);
            if( !error.empty() ){
                tmpError.push_back(error[i]);
            }
        }
    }

    CSpectrumSpectralAxis *_SpectralAxis = new CSpectrumSpectralAxis(tmpWave.size(), otherSpectral.IsInLogScale());
    CSpectrumFluxAxis *_FluxAxis = new CSpectrumFluxAxis(tmpFlux.size());

    for(Int32 i=0; i<tmpFlux.size(); i++){
        (*_SpectralAxis)[i] = tmpWave[i];
        (*_FluxAxis)[i] = tmpFlux[i];
        if( !error.empty() ){
            (*_FluxAxis).GetError()[i] = tmpError[i];
        }
    }

    m_SpectralAxis = *_SpectralAxis;
    m_FluxAxis = *_FluxAxis;

    m_estimationMethod = other.GetContinuumEstimationMethod();
    m_dfBinPath = other.GetWaveletsDFBinPath();
    m_medianWindowSize = other.GetMedianWinsize();
    m_nbScales = other.GetDecompScales();
}

CSpectrum::CSpectrum(CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis) :
  m_SpectralAxis(spectralAxis),
  m_FluxAxis(fluxAxis)
{
    m_estimationMethod = "";
    m_medianWindowSize = -1;
    m_nbScales = -1;
    m_dfBinPath = "";
}

CSpectrum::~CSpectrum()
{

}

CSpectrum& CSpectrum::operator=(const CSpectrum& other)
{
    m_SpectralAxis = other.GetSpectralAxis();
    m_FluxAxis = other.GetFluxAxis();

    m_estimationMethod = other.GetContinuumEstimationMethod();
    m_dfBinPath = other.GetWaveletsDFBinPath();
    m_medianWindowSize = other.GetMedianWinsize();
    m_nbScales = other.GetDecompScales();

    return *this;
}


Bool CSpectrum::RemoveContinuum(  CContinuum& remover )
{
    CSpectrumFluxAxis fluxAxisWithoutContinuum;

    Bool ret = remover.RemoveContinuum( *this, fluxAxisWithoutContinuum );

    m_FluxAxis = fluxAxisWithoutContinuum;

    return ret;
}

/**
 * Invert the flux axis
 */
Bool CSpectrum::InvertFlux()
{
    m_FluxAxis.Invert();
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

bool CSpectrum::GetMeanAndStdFluxInRange(TFloat64Range wlRange,  Float64& mean, Float64 &std) const
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
    const TFloat64List& error = m_FluxAxis.GetError();
    Float64 _Mean = 0.0;
    Float64 _SDev = 0.0;
    m_FluxAxis.ComputeMeanAndSDev( mask,_Mean ,_SDev, error);

    mean = _Mean;
    std = _SDev;
    return true;
}

bool CSpectrum::GetLinearRegInRange(TFloat64Range wlRange,  Float64& a, Float64 &b) const
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

    const TFloat64List& error = m_FluxAxis.GetError();

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
        y[k] = m_FluxAxis[ik];
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

void CSpectrum::SetName( const char* name )
{
    m_Name = name;
}

const Bool CSpectrum::checkFlux( Float64 flux, Int32 index ) const
{
    Log.LogDebug("    CSpectrum::checkFlux - Found flux value (=%e) at index=%d", flux, index);
    Bool invalidValue = false;
    if( std::isnan(flux) ){
        invalidValue = true;
        Log.LogDebug("    CSpectrum::checkFlux - Found nan flux value (=%e) at index=%d", flux, index);
    }
    if( std::isinf(flux) ){
        invalidValue = true;
        Log.LogDebug("    CSpectrum::checkFlux - Found inf flux value (=%e) at index=%d", flux, index);
    }
    if( flux != flux ){
        invalidValue = true;
        Log.LogDebug("    CSpectrum::checkFlux - Found invalid flux value (=%e) at index=%d", flux, index);
    }
    return invalidValue;
}

const Bool CSpectrum::checkNoise( Float64 error, Int32 index ) const
{
    Log.LogDebug("    CSpectrum::checkNoise - Found noise value (=%e) at index=%d", error, index);
    Bool validValue=true;
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

const Bool CSpectrum::IsFluxValid( Float64 LambdaMin,  Float64 LambdaMax ) const
{
    Bool allzero=true;
    Bool atleastOneInvalidValue = false;

    const Float64* flux = m_FluxAxis.GetSamples();
    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    Log.LogDetail( "CSpectrum::IsFluxValid - checking on the configured lambdarange = (%f, %f)", LambdaMin, LambdaMax );
    Log.LogDetail( "CSpectrum::IsFluxValid - checking on the true observed spectral axis lambdarange = (%f, %f)", m_SpectralAxis[iMin], m_SpectralAxis[iMax] );
    for(Int32 i=iMin; i<iMax; i++){
        //check flux
        atleastOneInvalidValue = checkFlux(flux[i], i);

        //all zero check
        if( flux[i] != 0.0 ){
            allzero = false;
            //Log.LogDebug("    CSpectrum::IsFluxValid - Found non zero and valid flux value (=%e) at index=%d", i, flux[i]);
        }
    }
    Bool valid = !allzero && !atleastOneInvalidValue;
    return valid;
}

const Bool CSpectrum::IsNoiseValid( Float64 LambdaMin,  Float64 LambdaMax ) const
{
    Bool valid=true;
    Int32 nInvalid = 0;

    const TFloat64List& error = m_FluxAxis.GetError();
    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    Log.LogDebug("    CSpectrum::IsNoiseValid - wl=%f iMin=%d error[iMin]=%e", LambdaMin, iMin, error[iMin]);
    Log.LogDebug("    CSpectrum::IsNoiseValid - wl=%f iMax=%d error[iMax]=%e", LambdaMax, iMax, error[iMax]);
    for(Int32 i=iMin; i<iMax; i++){
        //check noise
        Bool validSample = checkNoise(error[i], i);

        if(!validSample){
            valid=false;
            nInvalid++;
        }
    }
    if(nInvalid>0)
    {
        Log.LogDetail("    CSpectrum::IsNoiseValid - Found %d invalid noise samples", nInvalid);
    }
    return valid;
}

Bool CSpectrum::correctSpectrum( Float64 LambdaMin,  Float64 LambdaMax, Float64 coeffCorr )
{
    Bool corrected=false;
    Int32 nCorrected = 0;

    TFloat64List& error = m_FluxAxis.GetError();
    Float64* flux = m_FluxAxis.GetSamples();

    Int32 iMin = m_SpectralAxis.GetIndexAtWaveLength(LambdaMin);
    Int32 iMax = m_SpectralAxis.GetIndexAtWaveLength(LambdaMax);
    //Log.LogDebug("    CSpectrum::correctSpectrum - debug - iMin=%d and wmin=%f, iMax=%d and wmax=%f", iMin, m_SpectralAxis[iMin], iMax, m_SpectralAxis[iMax]);

    Float64 maxNoise = -DBL_MAX;
    Float64 minFlux = DBL_MAX;
    for(Int32 i=iMin; i<=iMax; i++){
        //Log.LogDebug("    CSpectrum::correctSpectrum - debug - RAW sample for maxFlux/minNoise. Found err(=%f) and flux=%f", error[i], flux[i]);
        //check noise
        if( error[i] < DBL_MIN ){
            continue;
        }
        if( std::isnan(error[i]) ){
            continue;
        }
        if( std::isinf(error[i]) ){
            continue;
        }
        if( error[i] != error[i] ){
            continue;
        }

        //check flux
        if( std::isnan(flux[i]) ){
            continue;
        }
        if( std::isinf(flux[i]) ){
            continue;
        }
        if( flux[i] != flux[i] ){
            continue;
        }
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
        return false;
    }

    for(Int32 i=iMin; i<=iMax; i++){
        Bool validSample=true;

        //check noise
        validSample = checkNoise(error[i], i);

        //check flux
        validSample = !checkFlux(flux[i], i);

        if(!validSample){
            error[i]=maxNoise*coeffCorr;
            flux[i]=minFlux/coeffCorr;
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

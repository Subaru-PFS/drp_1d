#ifndef _REDSHIFT_SPECTRUM_SPECTRUM_
#define _REDSHIFT_SPECTRUM_SPECTRUM_

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/continuum/continuum.h>

#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CSpectrum
{

public:

    enum EFLags
    {

    };

    CSpectrum();
    CSpectrum(const CSpectrum& other, TFloat64List mask);
    CSpectrum(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis);
    ~CSpectrum();

    CSpectrum& operator=(const CSpectrum& other);

    void  SetName( const char* name );

    Bool InvertFlux();

    const CSpectrumSpectralAxis&    GetSpectralAxis() const;
    const CSpectrumFluxAxis&        GetFluxAxis() const;

    const std::string               GetName() const;

    CSpectrumFluxAxis&              GetFluxAxis();
    CSpectrumSpectralAxis&          GetSpectralAxis();

    UInt32                          GetSampleCount() const;
    Float64                         GetResolution() const;
    Float64                         GetMeanResolution() const;
    TLambdaRange                    GetLambdaRange() const;

    bool                            GetMeanAndStdFluxInRange( TFloat64Range wlRange, Float64& mean, Float64& std ) const;
    bool                            GetLinearRegInRange( TFloat64Range wlRange,  Float64& a, Float64& b) const;

    Bool                            ConvertToLogScale();
    Bool                            ConvertToLinearScale();

    Bool                            RemoveContinuum( CContinuum& remover );
    const Bool                      checkFlux(Float64 flux, Int32 index) const;
    const Bool                      checkNoise(Float64 error, Int32 index) const;
    const Bool                      IsFluxValid(Float64 LambdaMin, Float64 LambdaMax) const;
    const Bool                      IsNoiseValid(Float64 LambdaMin, Float64 LambdaMax) const;
    Bool                            correctSpectrum(Float64 LambdaMin, Float64 LambdaMax, Float64 coeffCorr=10.0);

    const std::string&       	    GetFullPath() const;
    const Int32                     GetDecompScales() const;
    const Float64                   GetMedianWinsize() const;
    const std::string               GetContinuumEstimationMethod() const;
    const std::string               GetWaveletsDFBinPath() const;

    void 			    SetFullPath(const char* nameP);
    void 			    SetDecompScales(Int32 decompScales);
    void 			    SetMedianWinsize(Float64 winsize);
    void                SetContinuumEstimationMethod(std::string method);
    void                SetWaveletsDFBinPath(std::string binPath);

    void                LoadSpectrum(const char* spectrumFilePath, const char* noiseFilePath);

    Bool                SetTemplateBuffer( CSpectrum spectrum);
    /*Bool              Rebin( const TFloat64Range& range, const CSpectrumSpectralAxis& targetSpectralAxis,
                               CSpectrum& rebinedSpectrum, CMask& rebinedMask, const std::string opt_interp ); //linear always
                               */
    Bool                Rebin2( const TFloat64Range& range, const CSpectrumSpectralAxis& targetSpectralAxis,
                                CSpectrumFluxAxis& rebinedFluxAxis, CSpectrumSpectralAxis& rebinedSpectralAxis,/*CSpectrum& rebinedSpectrum,*/ CMask& rebinedMask, const std::string opt_interp = "lin", Float64 sourcez = -1 ); 
                               //sourcez only used with precomputefinegrid  
protected:
    CSpectrumSpectralAxis           m_SpectralAxis;
    CSpectrumFluxAxis               m_FluxAxis;
    TFloat64List                    m_pfgTplBuffer;
private:

    std::string                     m_Name;
    std::string                     m_FullPath;
    Int32                           m_nbScales;
    Float64                         m_medianWindowSize;
    std::string                     m_estimationMethod;
    std::string                     m_dfBinPath;
};

inline
UInt32 CSpectrum::GetSampleCount() const
{
    return m_SpectralAxis.GetSamplesCount();
}

inline
const CSpectrumSpectralAxis& CSpectrum::GetSpectralAxis() const
{
    return m_SpectralAxis;
}

inline
const CSpectrumFluxAxis& CSpectrum::GetFluxAxis() const
{
    return m_FluxAxis;
}

inline
CSpectrumSpectralAxis& CSpectrum::GetSpectralAxis()
{
    return m_SpectralAxis;
}

inline
CSpectrumFluxAxis& CSpectrum::GetFluxAxis()
{
    return m_FluxAxis;
}


}

#endif

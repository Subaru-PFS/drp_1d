#ifndef _REDSHIFT_SPECTRUM_SPECTRUM_
#define _REDSHIFT_SPECTRUM_SPECTRUM_

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/continuum/continuum.h>

#include <unordered_map>
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

    // FluxAxis components
    enum EType
    {
        nType_raw = 1,
        nType_continuumOnly = 2,
        nType_noContinuum = 3
    };

    CSpectrum();
    CSpectrum(const CSpectrum& other, TFloat64List mask);
    CSpectrum(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis);
    CSpectrum(const CSpectrum& other);
    ~CSpectrum();

    CSpectrum& operator=(const CSpectrum& other);

    void SetName(const char* name);
    void SetType(const EType type) const;

    const std::string&              GetName() const;
    const EType                     GetType() const;

    Bool InvertFlux();

    const CSpectrumSpectralAxis&    GetSpectralAxis() const;
    const CSpectrumFluxAxis&        GetFluxAxis() const;
    const CSpectrumFluxAxis&        GetRawFluxAxis() const;
    const CSpectrumFluxAxis&        GetContinuumFluxAxis() const;
    const CSpectrumFluxAxis&        GetWithoutContinuumFluxAxis() const;

    CSpectrumSpectralAxis&          GetSpectralAxis();
    CSpectrumFluxAxis&              GetFluxAxis();
    CSpectrumFluxAxis&              GetRawFluxAxis();
    CSpectrumFluxAxis&              GetContinuumFluxAxis();
    CSpectrumFluxAxis&              GetWithoutContinuumFluxAxis();

    UInt32                          GetSampleCount() const;
    Float64                         GetResolution() const;
    Float64                         GetMeanResolution() const;
    TLambdaRange                    GetLambdaRange() const;
    const std::string&              GetBaseline() const;

    bool                            GetMeanAndStdFluxInRange( TFloat64Range wlRange, Float64& mean, Float64& std ) const;
    bool                            GetLinearRegInRange( TFloat64Range wlRange,  Float64& a, Float64& b) const;

    Bool                            ConvertToLogScale();
    Bool                            ConvertToLinearScale();

    Bool                            RemoveContinuum( CContinuum& remover ) const;
    const Bool                      checkFlux(Float64 flux, Int32 index) const;
    const Bool                      checkNoise(Float64 error, Int32 index) const;
    const Bool                      IsFluxValid(Float64 LambdaMin, Float64 LambdaMax) const;
    const Bool                      IsNoiseValid(Float64 LambdaMin, Float64 LambdaMax) const;
    Bool                            correctSpectrum(Float64 LambdaMin, Float64 LambdaMax, Float64 coeffCorr=10.0);

    const std::string&       	    GetFullPath() const;
    const Int32                     GetDecompScales() const;
    const Float64                   GetMedianWinsize() const;
    const std::string&              GetContinuumEstimationMethod() const;
    const std::string&              GetWaveletsDFBinPath() const;

    void 			    SetFullPath(const char* nameP);
    void 			    SetDecompScales(Int32 decompScales);
    void 			    SetMedianWinsize(Float64 winsize);
    void                            SetContinuumEstimationMethod(std::string method);
    void                            SetContinuumEstimationMethod(const CSpectrumFluxAxis &ContinuumFluxAxis);
    void                            SetWaveletsDFBinPath(std::string binPath);

    void                            LoadSpectrum(const char* spectrumFilePath, const char* noiseFilePath);

    void                            ScaleFluxAxis(Float64 scale);

    void                            InitPrecomputeFineGrid() const;

    Bool                            Rebin( const TFloat64Range& range, const CSpectrumSpectralAxis& targetSpectralAxis,
                                           CSpectrum& rebinedSpectrum, CMask& rebinedMask, const std::string opt_interp = "lin",
                                           const std::string opt_error_interp="no" ) const;

protected:

    CSpectrumSpectralAxis           m_SpectralAxis;

    void                            EstimateContinuum() const;
    void                            ResetContinuum() const;
    Bool                            RebinFineGrid() const;

    const Float64                   m_dLambdaFineGrid = 0.1; //oversampling step for fine grid
                                                             //check if enough to be private
    mutable TFloat64List            m_pfgFlux;
    mutable Bool                    m_FineGridInterpolated = false;

    std::string                     m_Name;
    std::string                     m_FullPath;

    // Continuum removal parameters
    Int32                           m_nbScales;
    Float64                         m_medianWindowSize;
    std::string                     m_estimationMethod;
    std::string                     m_dfBinPath;

    mutable EType                   m_spcType = nType_raw;
    CSpectrumFluxAxis               m_RawFluxAxis;
    mutable CSpectrumFluxAxis       m_ContinuumFluxAxis;
    mutable CSpectrumFluxAxis       m_WithoutContinuumFluxAxis;

    // Flag
    mutable bool                    alreadyRemoved = false;

    // Map method2baseline
    const std::unordered_map<std::string, std::string> m_method2baseline = {
        {"IrregularSamplingMedian", "baselineISMedian"},
        {"Median",                  "baselineMedian"},
        {"waveletsDF",              "baselineDF"},
        {"raw",                     "baselineRAW"},
        {"zero",                    "baselineZERO"},
        {"manual",                  "baselineMANUAL"}
    };

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
    switch(m_spcType){
        case nType_raw :
                return GetRawFluxAxis();
                break;
        case nType_continuumOnly :
                return GetContinuumFluxAxis();
                break;
        case nType_noContinuum :
                return GetWithoutContinuumFluxAxis();
                break;
        default :
                return GetRawFluxAxis();
    }
}

inline
const CSpectrumFluxAxis& CSpectrum::GetRawFluxAxis() const
{
    return m_RawFluxAxis;
}

inline
const CSpectrumFluxAxis& CSpectrum::GetContinuumFluxAxis() const
{
    if( !alreadyRemoved ) {
        EstimateContinuum();
    }
    return m_ContinuumFluxAxis;
}

inline
const CSpectrumFluxAxis& CSpectrum::GetWithoutContinuumFluxAxis() const
{
    if( !alreadyRemoved ) {
        EstimateContinuum();
    }
    return m_WithoutContinuumFluxAxis;
}

inline
CSpectrumSpectralAxis& CSpectrum::GetSpectralAxis()
{
    return m_SpectralAxis;
}

inline
CSpectrumFluxAxis& CSpectrum::GetFluxAxis()
{
    switch(m_spcType){
        case nType_raw :
                return GetRawFluxAxis();
                break;
        case nType_continuumOnly :
                return GetContinuumFluxAxis();
                break;
        case nType_noContinuum :
                return GetWithoutContinuumFluxAxis();
                break;
        default :
                return GetRawFluxAxis();
    }
}

inline
CSpectrumFluxAxis& CSpectrum::GetRawFluxAxis()
{
    return m_RawFluxAxis;
}

inline
CSpectrumFluxAxis& CSpectrum::GetContinuumFluxAxis()
{
    if( !alreadyRemoved ) {
        EstimateContinuum();
    }
    return m_ContinuumFluxAxis;
}

inline
CSpectrumFluxAxis& CSpectrum::GetWithoutContinuumFluxAxis()
{
    if( !alreadyRemoved ) {
        EstimateContinuum();
    }
    return m_WithoutContinuumFluxAxis;
}

}

#endif

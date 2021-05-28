#ifndef _REDSHIFT_SPECTRUM_SPECTRUM_
#define _REDSHIFT_SPECTRUM_SPECTRUM_

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/spectrum/LSF.h>
#include <RedshiftLibrary/spectrum/LSFConstant.h>
#include <RedshiftLibrary/continuum/continuum.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
#include <unordered_map>
#include <stdexcept>
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
    CSpectrum(const std::string& name);
    CSpectrum(const CSpectrum& other, TFloat64List mask);
    CSpectrum(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis);
    CSpectrum(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, const std::shared_ptr<CLSF>& lsf);
    CSpectrum(const CSpectrum& other);
    ~CSpectrum();

    CSpectrum& operator=(const CSpectrum& other);

    void InitSpectrum(CParameterStore& parameterStore);
    void SetName(const std::string name);
    void SetType(const EType type) const;

    const std::string&              GetName() const;
    const EType                     GetType() const;

    Bool InvertFlux();

    const CSpectrumSpectralAxis&    GetSpectralAxis() const;
    const CSpectrumFluxAxis&        GetFluxAxis() const;
    const CSpectrumFluxAxis&        GetRawFluxAxis() const;
    const CSpectrumFluxAxis&        GetContinuumFluxAxis() const;
    const CSpectrumFluxAxis&        GetWithoutContinuumFluxAxis() const;
    const CSpectrumNoiseAxis&       GetErrorAxis() const;
    std::shared_ptr<const CLSF>     GetLSF() const;

    void                            SetSpectralAxis(const CSpectrumSpectralAxis & spectralaxis);
    void                            SetSpectralAxis(CSpectrumSpectralAxis && spectralaxis);
    void                            SetFluxAxis(const CSpectrumFluxAxis & fluxaxis);
    void                            SetFluxAxis(CSpectrumFluxAxis && fluxaxis);
    void                            SetSpectralAndFluxAxes(CSpectrumSpectralAxis spcaxis, CSpectrumFluxAxis fluxaxis);
    void                            SetErrorAxis(const CSpectrumNoiseAxis & noiseaxis);
    void                            SetErrorAxis(CSpectrumNoiseAxis && noiseaxis);

    std::shared_ptr<CLSF>           GetLSF();
    void                            SetLSF(const std::shared_ptr<CLSF>& lsf);

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

    void 			                SetFullPath(const char* nameP);
    void 			                SetDecompScales(Int32 decompScales);
    void 			                SetMedianWinsize(Float64 winsize);
    void                            SetContinuumEstimationMethod(std::string method);
    void                            SetContinuumEstimationMethod(const CSpectrumFluxAxis &ContinuumFluxAxis);
    void                            SetWaveletsDFBinPath(std::string binPath);

    void                            ScaleFluxAxis(Float64 scale);

    Bool                            Rebin( const TFloat64Range& range, const CSpectrumSpectralAxis& targetSpectralAxis,
                                           CSpectrum& rebinedSpectrum, CMask& rebinedMask, const std::string opt_interp = "lin",
                                           const std::string opt_error_interp="no" ) const;
    Int32                           extractFrom(const CSpectrum& other, Int32 startIdx, Int32 endIdx);

protected:

    // protected mutable getters
    CSpectrumFluxAxis&              GetFluxAxis_(); 
    CSpectrumFluxAxis&              GetRawFluxAxis_();
    CSpectrumFluxAxis&              GetContinuumFluxAxis_();
    CSpectrumFluxAxis&              GetWithoutContinuumFluxAxis_();

    CSpectrumSpectralAxis           m_SpectralAxis;
    std::shared_ptr<CLSF>           m_LSF;

    void                            EstimateContinuum() const;
    void                            ResetContinuum() const;
    Bool                            RebinFineGrid() const;
    void                            ClearFineGrid() const;


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
CSpectrumFluxAxis& CSpectrum::GetFluxAxis_() 
{
    switch(m_spcType){
        case nType_raw :
                return GetRawFluxAxis_();
                break;
        case nType_continuumOnly :
                return GetContinuumFluxAxis_();
                break;
        case nType_noContinuum :
                return GetWithoutContinuumFluxAxis_();
                break;
        default :
                return GetRawFluxAxis_();
    }
}

inline
const CSpectrumFluxAxis& CSpectrum::GetRawFluxAxis() const
{
    return m_RawFluxAxis;
}

inline
CSpectrumFluxAxis& CSpectrum::GetRawFluxAxis_()
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
CSpectrumFluxAxis& CSpectrum::GetContinuumFluxAxis_()
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
CSpectrumFluxAxis& CSpectrum::GetWithoutContinuumFluxAxis_()
{
    if( !alreadyRemoved ) {
        EstimateContinuum();
    }
    return m_WithoutContinuumFluxAxis;
}

inline
const CSpectrumNoiseAxis&  CSpectrum::GetErrorAxis() const
{
    return GetFluxAxis().GetError();
}

inline 
void CSpectrum::SetErrorAxis(const CSpectrumNoiseAxis & erroraxis)
{
    GetFluxAxis_().GetError() = erroraxis;
}

inline 
void CSpectrum::SetErrorAxis(CSpectrumNoiseAxis && erroraxis)
{
    GetFluxAxis_().GetError() = std::move(erroraxis);
}

inline
std::shared_ptr<const CLSF> CSpectrum::GetLSF() const
{
    return m_LSF;

}

inline
std::shared_ptr<CLSF> CSpectrum::GetLSF()
{
    return m_LSF;

}

inline 
void CSpectrum::SetLSF(const std::shared_ptr<CLSF>& lsf)
{
    m_LSF = lsf;
}

inline
Int32  CSpectrum::extractFrom(const CSpectrum& other, Int32 startIdx, Int32 endIdx)
{
    m_SpectralAxis.extractFrom(other.GetSpectralAxis(), startIdx, endIdx);
    m_RawFluxAxis.extractFrom(other.GetFluxAxis(), startIdx, endIdx);
    //probably we should do the save for all axis of the spectrum
    return 0;
}

}
#endif

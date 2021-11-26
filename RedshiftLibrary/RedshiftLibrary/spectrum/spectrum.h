// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_SPECTRUM_SPECTRUM_
#define _REDSHIFT_SPECTRUM_SPECTRUM_

#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/LSF.h"

#include "RedshiftLibrary/spectrum/LSFFactory.h"
//TODO: check if below are still required
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"

#include "RedshiftLibrary/continuum/continuum.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
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
    CSpectrum(const CSpectrum& other);
    CSpectrum(CSpectrum&& other);
    CSpectrum(const CSpectrum& other, const TFloat64List & mask);
    CSpectrum(CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis);
    CSpectrum(CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis, const std::shared_ptr<CLSF>& lsf);
    ~CSpectrum();

    CSpectrum& operator=(const CSpectrum& other);
    CSpectrum& operator=(CSpectrum&& other);

    void InitSpectrum(CParameterStore& parameterStore);
    void SetName(const std::string name);
    virtual void SetType(const EType type) const;
    virtual void SetType(const EType type) {const_cast<const CSpectrum *>(this)->SetType(type); };// non const version needed for the virtualization 
 
    const std::string&              GetName() const;
    const EType                     GetType() const;

    Bool InvertFlux();

    const CSpectrumSpectralAxis&    GetSpectralAxis() const;
    const CSpectrumFluxAxis&        GetFluxAxis() const;
    const CSpectrumFluxAxis&        GetRawFluxAxis() const;
    const CSpectrumFluxAxis&        GetContinuumFluxAxis() const;
    const CSpectrumFluxAxis&        GetWithoutContinuumFluxAxis() const;
    const CSpectrumNoiseAxis&       GetErrorAxis() const;
    const std::shared_ptr<const CLSF>     GetLSF() const;

    virtual void                    SetSpectralAxis(const CSpectrumSpectralAxis & spectralaxis);
    virtual void                    SetSpectralAxis(CSpectrumSpectralAxis && spectralaxis);
    virtual void                    SetFluxAxis(const CSpectrumFluxAxis & fluxaxis);
    virtual void                    SetFluxAxis(CSpectrumFluxAxis && fluxaxis);
    virtual void                    SetSpectralAndFluxAxes(CSpectrumSpectralAxis spcaxis, CSpectrumFluxAxis fluxaxis);
    void                            SetErrorAxis(const CSpectrumNoiseAxis & noiseaxis);
    void                            SetErrorAxis(CSpectrumNoiseAxis && noiseaxis);

    bool                            IsNoiseEmpty()const;
    bool                            IsEmpty() const;
    bool                            IsValid()const;
    void                            ValidateSpectrum(TFloat64Range lambdaRange, 
                                                    Bool enableInputSpcCorrect);
    void                            SetLSF(const std::shared_ptr<const CLSF>& lsf);

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
    const Float64                   GetMedianWinsize() const;
    const std::string&              GetContinuumEstimationMethod() const;

    void 			                SetFullPath(const char* nameP);
    void 			                SetMedianWinsize(Float64 winsize);
    void                            SetContinuumEstimationMethod(std::string method) const;
    void                            SetContinuumEstimationMethod(const CSpectrumFluxAxis &ContinuumFluxAxis);

    void                            ScaleFluxAxis(Float64 scale);

    Bool                            Rebin( const TFloat64Range& range, const CSpectrumSpectralAxis& targetSpectralAxis,
                                           CSpectrum& rebinedSpectrum, CMask& rebinedMask, const std::string opt_interp = "lin",
                                           const std::string opt_error_interp="no" ) const;
    CSpectrum                       extract(Int32 startIdx, Int32 endIdx)const;

protected:

    // protected mutable getters
    CSpectrumFluxAxis&              GetFluxAxis_(); 
    CSpectrumFluxAxis&              GetRawFluxAxis_();
    CSpectrumFluxAxis&              GetContinuumFluxAxis_();
    CSpectrumFluxAxis&              GetWithoutContinuumFluxAxis_();

    CSpectrumSpectralAxis           m_SpectralAxis;
    std::shared_ptr<const CLSF>           m_LSF;

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
    mutable Float64                         m_medianWindowSize;
    mutable std::string                     m_estimationMethod;

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
const std::shared_ptr<const CLSF> CSpectrum::GetLSF() const
{
    return m_LSF;
}

inline
bool CSpectrum::IsEmpty()const
{
    return m_SpectralAxis.isEmpty() || GetFluxAxis().isEmpty();
}

inline
bool CSpectrum::IsValid()const
{
    return m_SpectralAxis.GetSamplesCount() == GetFluxAxis().GetSamplesCount() && !IsEmpty();
}

inline
bool CSpectrum::IsNoiseEmpty()const
{
    return GetErrorAxis().isEmpty();
}
inline 
void CSpectrum::SetLSF(const std::shared_ptr<const CLSF>& lsf)
{
    m_LSF = lsf;
}

inline
CSpectrum  CSpectrum::extract(Int32 startIdx, Int32 endIdx) const
{
    CSpectrum copy = *this;
    copy.SetSpectralAndFluxAxes(copy.m_SpectralAxis.extract(startIdx, endIdx),
                                copy.GetFluxAxis().extract(startIdx, endIdx));
    return copy;
}

}
#endif

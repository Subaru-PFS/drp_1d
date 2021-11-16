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
#ifndef _REDSHIFT_SPECTRUM_SPECTRALAXIS_
#define _REDSHIFT_SPECTRUM_SPECTRALAXIS_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/axis.h"
#include "RedshiftLibrary/common/range.h"

#include <vector>

namespace NSEpic
{

class CMask;

/**
 * \ingroup Redshift
 */
class CSpectrumSpectralAxis : public CSpectrumAxis
{

public:

    enum EFlags
    {
        nFLags_LogScale = 1 << 0,
        nFLags_LogSampled = 1 << 1 //using second bit(shifting once): 0000 0010
    };

    enum EShiftDirection
    {
        nShiftForward =0,
        nShiftBackward
    };

    CSpectrumSpectralAxis();
    CSpectrumSpectralAxis( UInt32 n, Bool isLogScale =false);
    CSpectrumSpectralAxis( const TFloat64List & samples, Bool isLogScale=false, std::string AirVacuum="" );
    CSpectrumSpectralAxis( TFloat64List && samples, Bool isLogScale=false, std::string AirVacuum="" );
    CSpectrumSpectralAxis( const Float64* samples, UInt32 n, std::string AirVacuum="");
    CSpectrumSpectralAxis( const CSpectrumSpectralAxis& origin, Float64 redshift, EShiftDirection direction );

    Float64             GetResolution( Float64 atWavelength = -1.0 ) const;
    Float64             GetMeanResolution() const;

    void                ShiftByWaveLength(  const CSpectrumSpectralAxis& origin, Float64 wavelengthOffset, EShiftDirection direction );
    void                ShiftByWaveLength( Float64 wavelengthOffset, EShiftDirection direction );

    void                ApplyOffset(Float64 wavelengthOffset);

    Int32               GetIndexAtWaveLength( Float64 waveLength ) const;
    TInt32Range         GetIndexesAtWaveLengthRange( const TFloat64Range& waveLengthRange ) const;


    Bool                ConvertToLinearScale();
    Bool                ConvertToLogScale();
    Bool                IsInLogScale() const;
    Bool                IsInLinearScale() const;

    TLambdaRange        GetLambdaRange() const;
    Bool                ClampLambdaRange( const TFloat64Range& range, TFloat64Range& clampedRange ) const;
    void                GetMask( const TFloat64Range& range,  CMask& mask ) const;
    Float64             IntersectMaskAndComputeOverlapRate( const TFloat64Range& lambdaRange,  const CMask& omask ) const;
    void                SetLogScale();
    Bool                CheckLoglambdaSampling()const;
    Bool                IsLogSampled(Float64 logGridstep)const;
    Bool                IsLogSampled()const;
    Float64             GetlogGridStep() const;
    void                RecomputePreciseLoglambda();
    TFloat64List        GetSubSamplingMask(UInt32 ssratio) const;
    TFloat64List        GetSubSamplingMask(UInt32 ssratio, TFloat64Range lambdarange) const;
    TFloat64List        GetSubSamplingMask(UInt32 ssratio, const TInt32Range & ilbda) const;
    UInt32              GetLogSamplingIntegerRatio(Float64 logstep, Float64& modulo) const;
    
private:

    mutable UInt32      m_SpectralFlags = 0;
    mutable Bool        m_regularLogSamplingChecked=false;
    mutable Float64     m_regularLogSamplingStep = NAN; //sampling log step with which sampling was validated in CheckLoglambdaSampling 

};

}

#endif

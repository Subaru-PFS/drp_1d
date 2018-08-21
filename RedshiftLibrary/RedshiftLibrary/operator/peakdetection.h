#ifndef _REDSHIFT_OPERATOR_PEAKDETECTION_
#define _REDSHIFT_OPERATOR_PEAKDETECTION_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

#include <RedshiftLibrary/spectrum/spectralaxis.h>

namespace NSEpic
{

class CSpectrum;
class CSpectrumAxis;
class CPeakDetectionResult;

class CPeakDetection
{

public:

    CPeakDetection( Float64 windowSize = 250.0, Float64 cut = 5.0, UInt32 medianSmoothHalfWidth = 1, UInt32 enlargeRate = 2.0 , Float64 detectionnoiseoffset=0.0);
    ~CPeakDetection();

    std::shared_ptr<const CPeakDetectionResult> Compute( const CSpectrum& spectrum, const TLambdaRange& lambdaRange);

private:
    Float64 m_winsize;
    Float64 m_cut;
    UInt32 m_medianSmoothHalfWidth;
    UInt32 m_enlargeRate;
    Float64 m_detectionnoiseoffset;

    void FindPossiblePeaks(const CSpectrumAxis& smoothedFluxAxis, const CSpectrumSpectralAxis& spectralAxis, TInt32RangeList& peakList );
    void RedefineBorders( TInt32RangeList& peakList, const CSpectrumAxis& waves, const CSpectrumAxis& smoothFluxAxis, const CSpectrumAxis& fluxAxis );
    TInt32Range FindGaussianFitStartAndStop( Int32 i, const TInt32RangeList& peaksBorders, UInt32 enlargeRate, Int32 len );
    TInt32Range LimitGaussianFitStartAndStop(Int32 i, const TInt32RangeList& peaksBorders, Int32 len , const CSpectrum &spectrum);
    Float64 XMad( const Float64* x, Int32 n, Float64 median );

};

}

#endif

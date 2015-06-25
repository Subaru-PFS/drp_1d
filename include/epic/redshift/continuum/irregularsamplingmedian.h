#ifndef _REDSHIFT_CONTINUUM_IRREGULARSMAPLINGMEDIAN_
#define _REDSHIFT_CONTINUUM_IRREGULARSMAPLINGMEDIAN_

#include <epic/redshift/continuum/continuum.h>

namespace NSEpic
{

class CSpectrumFluxAxis;

class CContinuumIrregularSamplingMedian : public CContinuum
{

    DEFINE_MANAGED_OBJECT( CContinuumIrregularSamplingMedian )

public:

    CContinuumIrregularSamplingMedian();
    ~CContinuumIrregularSamplingMedian();

    Void SetMeanKernelWidth( Float32 width );
    Void SetMedianKernelWidth( Float32 width );
    Void SetMedianCycleCount( UInt32 count );

    Bool RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis );
    Bool ProcessRemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis, Float64 resolution );


private:

    Int32   MedianSmooth( const Float64 *y, Int32 n_points, Int32 n_range, Float64 *y_out );

    Int32   MeanSmooth( const Float64 *y, Int32 N, Int32 n, Float64 *y_out );

    Int32   OddMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64* y_out );
    Int32   EvenMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64* y_out );

    Int32   m_MeanSmoothAmplitude;
    Int32   m_MedianSmoothCycles;
    Int32   m_MedianSmoothAmplitude;
    Bool    m_Even;

};


}

#endif

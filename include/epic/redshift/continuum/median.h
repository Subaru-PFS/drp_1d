#ifndef _REDSHIFT_CONTINUUM_MEDIAN_
#define _REDSHIFT_CONTINUUM_MEDIAN_

#include <epic/redshift/continuum/continuum.h>

namespace __NS__
{

class CSpectrumFluxAxis;

class CContinuumMedian : public CContinuum
{

    DEFINE_MANAGED_OBJECT( CContinuumMedian )

public:

    CContinuumMedian();
    ~CContinuumMedian();

    Bool RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis );

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

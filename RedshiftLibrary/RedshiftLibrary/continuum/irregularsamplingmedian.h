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
#ifndef _REDSHIFT_CONTINUUM_IRREGULARSMAPLINGMEDIAN_
#define _REDSHIFT_CONTINUUM_IRREGULARSMAPLINGMEDIAN_

#include "RedshiftLibrary/continuum/continuum.h"

namespace NSEpic
{

class CSpectrumFluxAxis;

/** \ingroup Redshift
 * Algorithm for estimating the continuum by computing the 'medium' resolution and applying the median method to it.
 */
class CContinuumIrregularSamplingMedian : public CContinuum
{

public:

    CContinuumIrregularSamplingMedian();
    ~CContinuumIrregularSamplingMedian();

    void SetMeanKernelWidth( Float32 width );
    void SetMedianKernelWidth( Float32 width );
    void SetMedianCycleCount( UInt32 count );

    bool RemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis );
    bool ProcessRemoveContinuum( const CSpectrum& s, CSpectrumFluxAxis& noContinuumFluxAxis, Float64 resolution );


private:

    Int32   MedianSmooth( const Float64 *y, Int32 n_points, Int32 n_range, Float64 *y_out );

    Int32   MeanSmooth( const Float64 *y, Int32 N, Int32 n, Float64 *y_out );

    Int32   OddMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64 y_input_begin_val, Float64 y_input_end_val, Float64* y_out );
    Int32   EvenMirror( const Float64* y_input, Int32 N, Int32 Nreflex, Float64* y_out );

    Int32   m_MeanSmoothAmplitude;
    Int32   m_MedianSmoothCycles;
    Int32   m_MedianSmoothAmplitude;
    bool    m_Even;

};


}

#endif

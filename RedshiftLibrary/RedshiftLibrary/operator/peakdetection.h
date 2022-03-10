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
#ifndef _REDSHIFT_OPERATOR_PEAKDETECTION_
#define _REDSHIFT_OPERATOR_PEAKDETECTION_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include "RedshiftLibrary/spectrum/spectralaxis.h"

#include <iterator>

namespace test_peakdetection{ //boost_test_suite
    //all boost_auto_test_case that use private method
    class XMad_test;
    class FindGaussianFitStartAndStop_test;
    class FindPossiblePeaks_test;
    class RedefineBorders_test;
}
namespace NSEpic
{

class CSpectrum;
class CSpectrumAxis;
class CPeakDetectionResult;

class CPeakDetection
{

public:

    CPeakDetection( Float64 windowSize = 250.0, Float64 cut = 5.0, Int32 medianSmoothHalfWidth = 1, Int32 enlargeRate = 2.0 , Float64 detectionnoiseoffset=0.0);
    ~CPeakDetection();

    std::shared_ptr<const CPeakDetectionResult> Compute( const CSpectrum& spectrum);

private:
    friend class test_peakdetection::XMad_test;
    friend class test_peakdetection::FindGaussianFitStartAndStop_test;
    friend class test_peakdetection::FindPossiblePeaks_test;
    friend class test_peakdetection::RedefineBorders_test;

    Float64 m_winsize;
    Float64 m_cut;
    Int32 m_medianSmoothHalfWidth;
    Int32 m_enlargeRate;
    Float64 m_detectionnoiseoffset;

    void FindPossiblePeaks(const CSpectrumAxis& smoothedFluxAxis, const CSpectrumSpectralAxis& spectralAxis, TInt32RangeList& peakList );
    void RedefineBorders( TInt32RangeList& peakList, const CSpectrumAxis& waves, const CSpectrumAxis& smoothFluxAxis, const CSpectrumAxis& fluxAxis );
    TInt32Range FindGaussianFitStartAndStop( Int32 i, const TInt32RangeList& peaksBorders, Int32 enlargeRate, Int32 len );
    Float64 XMad( const TFloat64List::const_iterator &begin, const TFloat64List::const_iterator &end, Float64 median );
    Float64 XMad( const TFloat64List &x, Float64 median );

};

}

#endif

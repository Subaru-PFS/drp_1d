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
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/operator/peakdetection.h"
#include "RedshiftLibrary/operator/peakdetectionresult.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <numeric>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(test_peakdetection)

BOOST_AUTO_TEST_CASE(XMad_test) {
  CPeakDetection peakDetection;
  TFloat64List x;
  Float64 returnValue;

  CMedian<Float64> medianProcessor;
  Float64 med;

  // ODD RANGE
  x = {1., 2., 3., 2., 1.};

  med = medianProcessor.Find(x);
  BOOST_CHECK_CLOSE(med, 2., 1e-12);

  returnValue = peakDetection.XMad(x, med);
  BOOST_CHECK_CLOSE(returnValue, 1., 1e-12);

  // EVEN RANGE
  x = {1., 2., 3., 3., 2., 1.};

  med = medianProcessor.Find(x);
  BOOST_CHECK_CLOSE(med, 2., 1e-12);

  returnValue = peakDetection.XMad(x, med);
  BOOST_CHECK_CLOSE(returnValue, 1., 1e-12);

  Float64 returnValue_2 = peakDetection.XMad(x.begin(), x.end(), med);
  BOOST_CHECK_CLOSE(returnValue, returnValue_2, 1e-12);
}

BOOST_AUTO_TEST_CASE(FindGaussianFitStartAndStop_test) {
  Int32 n = 250;

  TInt32RangeList peak;
  peak.push_back(TInt32Range(-5, 35));
  peak.push_back(TInt32Range(30, 70));
  peak.push_back(TInt32Range(80, 120));
  peak.push_back(TInt32Range(160, 200));
  peak.push_back(TInt32Range(190, 230));

  CPeakDetection peakDetection;

  // len = 250
  // width = 41
  // enlargeRate = 1
  //
  // fitStart = max( 0, fitStart - (int)(enlargeRate*width) ); --> 0
  // fitStop = min( len, fitStop + (int)(enlargeRate*width) ); --> 77
  //
  // fitStop = min( peaksBorders[i+1].GetBegin(), fitStop ); --> 30
  TInt32Range range = peakDetection.FindGaussianFitStartAndStop(0, peak, 1, n);
  BOOST_CHECK_EQUAL(range, TInt32Range(0, 30));

  // len = 250
  // width = 41
  // enlargeRate = 1
  //
  // fitStart = max( 0, fitStart - (int)(enlargeRate*width) ); --> 0
  // fitStop = min( len, fitStop + (int)(enlargeRate*width) ); --> 112
  //
  // fitStart = max( peaksBorders[i-1].GetEnd(), fitStart ); --> 35
  // fitStop = min( peaksBorders[i+1].GetBegin(), fitStop ); --> 80
  range = peakDetection.FindGaussianFitStartAndStop(1, peak, 1, n);
  BOOST_CHECK_EQUAL(range, TInt32Range(35, 80));

  // len = 250
  // width = 41
  // enlargeRate = 1
  //
  // fitStart = max( 0, fitStart - (int)(enlargeRate*width) ); --> 39
  // fitStop = min( len, fitStop + (int)(enlargeRate*width) ); --> 162
  //
  // fitStart = max( peaksBorders[i-1].GetEnd(), fitStart ); --> 70
  // fitStop = min( peaksBorders[i+1].GetBegin(), fitStop ); --> 160
  range = peakDetection.FindGaussianFitStartAndStop(2, peak, 1, n);
  BOOST_CHECK_EQUAL(range, TInt32Range(70, 160));

  // len = 250
  // width = 41
  // enlargeRate = 1
  //
  // fitStart = max( 0, fitStart - (int)(enlargeRate*width) ); --> 119
  // fitStop = min( len, fitStop + (int)(enlargeRate*width) ); --> 242
  //
  // fitStart = max( peaksBorders[i-1].GetEnd(), fitStart ); --> 120
  // fitStop = min( peaksBorders[i+1].GetBegin(), fitStop ); --> 190
  range = peakDetection.FindGaussianFitStartAndStop(3, peak, 1, n);
  BOOST_CHECK_EQUAL(range, TInt32Range(120, 190));

  // len = 250
  // width = 41
  // enlargeRate = 1
  //
  // fitStart = max( 0, fitStart - (int)(enlargeRate*width) ); --> 149
  // fitStop = min( len, fitStop + (int)(enlargeRate*width) ); --> 250
  //
  // fitStart = max( peaksBorders[i-1].GetEnd(), fitStart ); --> 200
  range = peakDetection.FindGaussianFitStartAndStop(4, peak, 1, n);
  BOOST_CHECK_EQUAL(range, TInt32Range(200, 250));
}

BOOST_AUTO_TEST_CASE(FindPossiblePeaks_test) {
  CPeakDetection peakDetection;

  TFloat64List inc(10);
  std::iota(inc.begin(), inc.end(), 0.0);
  CSpectrumSpectralAxis spectralAxis(std::move(inc), false);

  CSpectrumFluxAxis modelfluxAxis(
      TFloat64List{1., 1., 1., 10., 20., 10., 1., 1., 1., 1.});
  CSpectrum spc(spectralAxis, modelfluxAxis);
  TInt32RangeList peakList;
  peakDetection.FindPossiblePeaks(modelfluxAxis, spectralAxis, peakList);

  BOOST_CHECK_EQUAL(peakList.size(), 1);
  BOOST_CHECK_EQUAL(peakList[0], TInt32Range(3, 5));

  // no peak detected
  Float64 *modelSamples = modelfluxAxis.GetSamples();
  for (Int32 k = 0; k < modelfluxAxis.GetSamplesCount(); k++) {
    modelSamples[k] = 1.;
  }
  spc.SetSpectralAndFluxAxes(spectralAxis, modelfluxAxis);
  peakDetection.FindPossiblePeaks(modelfluxAxis, spectralAxis, peakList);

  BOOST_CHECK_EQUAL(peakList.size(), 0);
}

BOOST_AUTO_TEST_CASE(RedefineBorders_test) {
  CPeakDetection peakDetection;

  TFloat64List inc(10);
  std::iota(inc.begin(), inc.end(), 0.0);
  CSpectrumSpectralAxis spectralAxis(std::move(inc), false);

  CSpectrumFluxAxis modelfluxAxis(
      TFloat64List{1., 1., 1., 10., 20., 10., 1., 1., 1., 1.});
  CSpectrum spc(spectralAxis, modelfluxAxis);

  TInt32RangeList peakList;
  peakDetection.FindPossiblePeaks(modelfluxAxis, spectralAxis, peakList);

  peakDetection.RedefineBorders(peakList, spectralAxis, modelfluxAxis,
                                modelfluxAxis);
  BOOST_CHECK_EQUAL(peakList.size(), 1);
  BOOST_CHECK_EQUAL(peakList[0].GetBegin(), 3);
  BOOST_CHECK_EQUAL(peakList[0].GetEnd(), 5);
}

BOOST_AUTO_TEST_CASE(Compute_test) {
  CPeakDetection peakDetection;

  TFloat64List inc(10);
  std::iota(inc.begin(), inc.end(), 0.0);
  CSpectrumSpectralAxis spectralAxis(std::move(inc), false);

  CSpectrumFluxAxis modelfluxAxis(
      TFloat64List{1., 1., 1., 10., 20., 10., 1., 1., 1., 1.});
  CSpectrum spc(spectralAxis, modelfluxAxis);

  std::shared_ptr<const CPeakDetectionResult> result;
  result = peakDetection.Compute(spc);
  BOOST_CHECK(result != NULL);
  BOOST_CHECK_EQUAL(result->PeakList[0], TInt32Range(3, 5));
  BOOST_CHECK_EQUAL(result->EnlargedPeakList[0], TInt32Range(0, 10));

  // no peak detected
  Float64 *modelSamples = modelfluxAxis.GetSamples();
  for (Int32 k = 0; k < modelfluxAxis.GetSamplesCount(); k++) {
    modelSamples[k] = 1.;
  }
  spc.SetSpectralAndFluxAxes(spectralAxis, modelfluxAxis);
  result = peakDetection.Compute(spc);

  BOOST_CHECK(result == NULL);
}

BOOST_AUTO_TEST_SUITE_END()
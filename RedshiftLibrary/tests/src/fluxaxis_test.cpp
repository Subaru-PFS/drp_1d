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
#include "RedshiftLibrary/spectrum/fluxaxis.h"

#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/mean.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/log/log.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <numeric>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Spectrum)

bool correctMessage(const GlobalException &ex) {
  BOOST_CHECK_EQUAL(
      ex.what(),
      std::string(
          "CSpectrum::Rebin: cannot interpolate outside input spectral range"));
  return true;
}
BOOST_AUTO_TEST_CASE(calcul) {

  //--------------------//
  // constructor

  CSpectrumFluxAxis object_FluxAxis;
  BOOST_CHECK(object_FluxAxis.GetSamplesCount() == 0);

  Int32 n = 10;
  CSpectrumFluxAxis object_FluxAxis2(n);
  BOOST_CHECK(object_FluxAxis2.GetSamplesCount() == n);

  Float64 Array1[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  CSpectrumFluxAxis object_FluxAxis3(Array1, 10);
  BOOST_CHECK(object_FluxAxis3.GetSamplesCount() == 10);
  BOOST_CHECK(accumulate(object_FluxAxis3.GetSamples(),
                         object_FluxAxis3.GetSamples() + 10, 0) == 55);
  BOOST_TEST_MESSAGE("index2:" << object_FluxAxis3.GetSamplesCount());

  Float64 Array2[10] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
  CSpectrumFluxAxis object_FluxAxis4(Array1, 10, Array2, 10);
  Int32 s = object_FluxAxis4.GetSamplesCount();
  BOOST_CHECK(object_FluxAxis4.GetSamplesCount() == 10);

  //----------//
  // test rbin

  TFloat64Range object_Range(1., 4.6);
  CSpectrumFluxAxis sourceFluxAxis(Array1, n);

  TFloat64List lbdaList(n);
  for (Int32 i = 0; i < n; i++) {
    lbdaList[i] = object_Range.GetBegin() + i * 0.4;
  }

  CSpectrumSpectralAxis sourceSpectralAxis(lbdaList);
  CSpectrumSpectralAxis targetSpectralAxis(lbdaList);

  CSpectrumFluxAxis rebinedFluxAxis(Array1, n);
  CSpectrumSpectralAxis rebinedSpectralAxis(lbdaList);
  CMask rebinedMask(n);

  // cas 1
  // std::shared_ptr<CSpectrum> object_CSpectrum = std::shared_ptr<CSpectrum>(
  // new CSpectrum( sourceSpectralAxis, sourceFluxAxis));
  CSpectrum object_CSpectrum(sourceSpectralAxis, sourceFluxAxis);
  CSpectrum rebinnedSpectrum;
  bool resultcas1 = object_CSpectrum.Rebin(object_Range, targetSpectralAxis,
                                           rebinnedSpectrum, rebinedMask);
  BOOST_CHECK(resultcas1 == true);
  BOOST_TEST_MESSAGE(
      "cas1:" << object_CSpectrum.Rebin(object_Range, targetSpectralAxis,
                                        rebinnedSpectrum, rebinedMask));

  // cas 3
  CSpectrumSpectralAxis sourceSpectralAxis3(lbdaList, 1);

  CSpectrum object_CSpectrum3(sourceSpectralAxis3, sourceFluxAxis);
  CSpectrum rebinnedSpectrum3;
  bool resultcas3 = object_CSpectrum3.Rebin(object_Range, targetSpectralAxis,
                                            rebinnedSpectrum3, rebinedMask);
  BOOST_CHECK(resultcas3 == false);
  BOOST_TEST_MESSAGE("cas3:" << resultcas3);

  // cas 4

  TFloat64Range object_Range4(3., 6.);
  TFloat64Range logIntersectedLambdaRange(log(object_Range4.GetBegin()),
                                          log(object_Range4.GetEnd()));
  TFloat64Range currentRange = logIntersectedLambdaRange;

  TFloat64List Array4 = {0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,
                         1.2,  1.3,  1.4,  1.45, 1.50, 1.55, 1.60,
                         1.65, 1.70, 1.75, 1.80, 1.85, 1.9};
  TFloat64List Array4b = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
  CSpectrumSpectralAxis targetSpectralAxis4(Array4, true);
  CSpectrumSpectralAxis sourceSpectralAxis4(Array4b, true);
  CMask rebinedMask2;
  CSpectrum object_CSpectrum4(sourceSpectralAxis4, sourceFluxAxis);
  CSpectrum rebinnedSpectrum4;
  /* BOOST_CHECK_EXCEPTION(
     object_CSpectrum4.Rebin(object_Range4,targetSpectralAxis4,rebinnedSpectrum4,(rebinedMask2)),
                         GlobalException, correctMessage);
 */
  //----------//
  // test rbin2

  // cas 1
  Float64 source = 1;
  const std::string opt_interp = "lin";

  CSpectrum object_CSpectrum5(sourceSpectralAxis, sourceFluxAxis);
  CMask rebinedMask3;
  CSpectrum rebinnedSpectrum5;
  bool resultRebincas1 =
      object_CSpectrum5.Rebin(object_Range, targetSpectralAxis,
                              rebinnedSpectrum5, (rebinedMask3), opt_interp);
  BOOST_CHECK(resultRebincas1 == true);
  BOOST_TEST_MESSAGE("Rebin cas1:" << resultRebincas1);

  // cas 3
  CMask rebinedMask5;
  CSpectrum object_CSpectrum7(sourceSpectralAxis3, sourceFluxAxis);
  CSpectrum rebinnedSpectrum7;
  bool resultRebincas3 =
      object_CSpectrum7.Rebin(object_Range, targetSpectralAxis,
                              rebinnedSpectrum7, (rebinedMask5), opt_interp);
  BOOST_CHECK(resultRebincas3 == false);
  BOOST_TEST_MESSAGE("Rebin cas3:" << resultRebincas3);

  // cas 4
  TFloat64List Array5b = {0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,
                          1.2,  1.3,  1.4,  1.45, 1.50, 1.55, 1.60,
                          1.65, 1.70, 1.75, 1.80, 1.85, 1.9};
  TFloat64List Array5 = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
  CSpectrumSpectralAxis targetSpectralAxis5(Array5, true);
  CSpectrumSpectralAxis sourceSpectralAxis5(Array5b, true);
  CSpectrumFluxAxis sourceFluxAxis5(Array1, n + n);

  CSpectrum object_CSpectrum8(sourceSpectralAxis5, sourceFluxAxis5);
  CMask rebinedMask6;
  CSpectrum rebinnedSpectrum8;
  bool resultRebincas4 =
      object_CSpectrum8.Rebin(object_Range4, targetSpectralAxis5,
                              rebinnedSpectrum8, (rebinedMask6), opt_interp);
  BOOST_CHECK(resultRebincas4 == true);
  BOOST_TEST_MESSAGE("Rebin cas4:" << resultRebincas4);

  // cas 5
  const std::string opt_interp2 = "precomputedfinegrid";

  CSpectrum object_CSpectrum9(sourceSpectralAxis5, sourceFluxAxis5);

  CMask rebinedMask7;
  CSpectrum rebinnedSpectrum9;
  bool resultRebincas5 =
      object_CSpectrum9.Rebin(object_Range4, targetSpectralAxis5,
                              rebinnedSpectrum9, rebinedMask7, opt_interp2);
  BOOST_CHECK(resultRebincas5 == true);
  BOOST_TEST_MESSAGE("Rebin cas5:" << resultRebincas5);

  CMask rebinedMask8;
  CSpectrum object_CSpectrum10(sourceSpectralAxis5, sourceFluxAxis5);
  CSpectrum object_CSpectrum_null;
  CSpectrum rebinnedSpectrum10;
  bool resultRebincas5bis =
      object_CSpectrum10.Rebin(object_Range4, targetSpectralAxis4,
                               rebinnedSpectrum10, rebinedMask8, opt_interp2);
  BOOST_CHECK(resultRebincas5bis == true);
  BOOST_TEST_MESSAGE("Rebin cas5bis:" << resultRebincas5bis);
  // cas 6
  const std::string opt_interp3 = "spline";

  CMask rebinedMask9;
  CSpectrum object_CSpectrum11(sourceSpectralAxis5, sourceFluxAxis5);
  CSpectrum rebinnedSpectrum11;
  bool resultRebincas6 =
      object_CSpectrum11.Rebin(object_Range4, targetSpectralAxis4,
                               rebinnedSpectrum11, rebinedMask9, opt_interp3);
  BOOST_CHECK(resultRebincas6 == true);
  BOOST_TEST_MESSAGE("Rebin cas6:" << resultRebincas6);

  // cas 7
  const std::string opt_interp4 = "ngp";
  CMask rebinedMask10;
  CSpectrum object_CSpectrum12(sourceSpectralAxis5, sourceFluxAxis5);
  CSpectrum rebinnedSpectrum12;
  bool resultRebincas7 = object_CSpectrum12.Rebin(
      object_Range4, targetSpectralAxis4, rebinnedSpectrum12, (rebinedMask10),
      opt_interp4);
  BOOST_CHECK(resultRebincas7 == true);
  BOOST_TEST_MESSAGE("Rebin cas7:" << resultRebincas7);

  //---------//
  // test ApplyMedianSmooth
  CSpectrumFluxAxis object_CSpectrumFluxAxis(Array1, n);

  bool resultApplyMedianSmoothcas1 =
      object_CSpectrumFluxAxis.ApplyMedianSmooth(0);
  BOOST_CHECK(resultApplyMedianSmoothcas1 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMedianSmooth cas1:" << resultApplyMedianSmoothcas1);

  BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis.GetSamplesCount"
                     << object_CSpectrumFluxAxis.GetSamplesCount());

  bool resultApplyMedianSmoothcas2 =
      object_CSpectrumFluxAxis.ApplyMedianSmooth(15);
  BOOST_CHECK(resultApplyMedianSmoothcas2 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMedianSmooth cas2:" << resultApplyMedianSmoothcas2);

  bool resultApplyMedianSmoothcas3 =
      object_CSpectrumFluxAxis.ApplyMedianSmooth(5);
  BOOST_CHECK(resultApplyMedianSmoothcas3 == true);
  BOOST_TEST_MESSAGE(
      "resultApplyMedianSmooth cas3:" << resultApplyMedianSmoothcas3);

  //--------------------//
  // ApplyMeanSmooth

  bool resultApplyMeanSmoothcas1 = object_CSpectrumFluxAxis.ApplyMeanSmooth(0);
  BOOST_CHECK(resultApplyMeanSmoothcas1 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMeanSmooth cas1:" << resultApplyMeanSmoothcas1);

  BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis.GetSamplesCount"
                     << object_CSpectrumFluxAxis.GetSamplesCount());

  bool resultApplyMeanSmoothcas2 = object_CSpectrumFluxAxis.ApplyMeanSmooth(15);
  BOOST_CHECK(resultApplyMeanSmoothcas2 == false);
  BOOST_TEST_MESSAGE(
      "resultApplyMeanSmooth cas2:" << resultApplyMeanSmoothcas2);

  bool resultApplyMeanSmoothcas3 = object_CSpectrumFluxAxis.ApplyMeanSmooth(5);
  BOOST_CHECK(resultApplyMeanSmoothcas3 == true);
  BOOST_TEST_MESSAGE(
      "resultApplyMeanSmooth cas3:" << resultApplyMeanSmoothcas3);

  //--------------------//
  // test ComputeRMSDiff

  Float64 ArrayA[] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  CSpectrumFluxAxis object_FluxAxisA(ArrayA, 10);

  Float64 ArrayB[] = {2., 4., 6., 8., 10., 12., 14., 16., 18., 20.};
  CSpectrumFluxAxis object_FluxAxisB(ArrayB, 10);

  Float64 resultComputeRMSDiff =
      object_FluxAxisA.ComputeRMSDiff(object_FluxAxisB);
  BOOST_TEST_MESSAGE(
      "object_FluxAxisA.GetSamples()[1]=" << object_FluxAxisA.GetSamples()[1]);
  BOOST_TEST_MESSAGE(
      "object_FluxAxisB.GetSamples()[1]=" << object_FluxAxisB.GetSamples()[1]);

  Float64 er2 = 0.f;
  Float64 er = 0.f;

  Float64 weight = (Float64)n;
  for (int j = 0; j < 10; j++) {
    er2 +=
        (object_FluxAxisA.GetSamples()[j] - object_FluxAxisB.GetSamples()[j]) *
        (object_FluxAxisA.GetSamples()[j] - object_FluxAxisB.GetSamples()[j]) /
        weight;
  }
  er = sqrt(er2);

  BOOST_CHECK_CLOSE(resultComputeRMSDiff, er, 1.e-12);

  //--------------------//
  // test Subtract

  bool resultSubtract = object_FluxAxisA.Subtract(object_FluxAxisB);

  int indice = 0;
  bool sub;

  for (Int32 i = 0; i < 10; i++) {
    object_FluxAxisA.GetSamples()[i] =
        object_FluxAxisA.GetSamples()[i] - object_FluxAxisB.GetSamples()[i];
    indice++;
  }

  if (indice == 10)
    sub = true;
  else
    sub = false;

  BOOST_CHECK(resultSubtract == sub);
  BOOST_TEST_MESSAGE("resultSubtract=" << resultSubtract << ", sub=" << sub);

  //--------------------//
  // test Invert

  bool resultInvert = object_FluxAxisA.Invert();

  int indice2 = 0;
  bool inv;

  for (Int32 i = 0; i < 10; i++) {
    object_FluxAxisA.GetSamples()[i] = -object_FluxAxisA.GetSamples()[i];

    indice2++;
  }

  if (indice2 == 10) {

    inv = true;
  } else {

    inv = false;
  }

  BOOST_CHECK(resultInvert == inv);
  BOOST_TEST_MESSAGE("resultInvert=" << resultInvert << ", inv=" << inv);

  //--------------------//
  // test ComputeMeanAndSDevWithoutError

  CMask Mask(10);

  Float64 mean = 1.;
  Float64 sdev = 1.;
  const CSpectrumNoiseAxis error(1);
  const CSpectrumNoiseAxis empty_error;

  object_FluxAxisA.GetError() = error;
  bool resultComputeMeanAndSDev_cas1 =
      object_FluxAxisA.ComputeMeanAndSDev(Mask, mean, sdev);
  BOOST_CHECK(resultComputeMeanAndSDev_cas1 == false);

  for (int i = 0; i < 10; i++) {
    Mask[i] = 0;
  }

  object_FluxAxisA.GetError() = empty_error;
  bool resultComputeMeanAndSDev_cas2 =
      object_FluxAxisA.ComputeMeanAndSDev(Mask, mean, sdev);
  BOOST_CHECK(resultComputeMeanAndSDev_cas2 == false);
}

BOOST_AUTO_TEST_CASE(RebinVarianceWeighted) {
  // test RebinVarianceWeighted
  TFloat64List lambdas = {1000, 2000, 3000, 4000, 5000,
                          6000, 7000, 8000, 9000, 10000};
  CSpectrumFluxAxis sourceFluxAxis(10);
  TFloat64List Array = {1000, 2500, 2900, 3890, 4690,
                        5500, 6800, 7001, 8033, 10000};
  CSpectrumSpectralAxis sourceSpectralAxis(Array, false);
  CSpectrumSpectralAxis bogus_sourceSpectralAxis(4, false);

  CSpectrumSpectralAxis targetSpectralAxis(lambdas);
  TFloat64Range currentRange = targetSpectralAxis.GetLambdaRange();

  const std::string opt_interp = "lin";
  std::string errorRebinMethod = "rebinVariance";

  CSpectrum object_CSpectrum13(sourceSpectralAxis, sourceFluxAxis);
  CMask rebinedMask13;
  CSpectrum rebinnedSpectrum13;
  bool resultRebinvaria1 = object_CSpectrum13.Rebin(
      currentRange, targetSpectralAxis, rebinnedSpectrum13, rebinedMask13,
      opt_interp, errorRebinMethod);
  BOOST_CHECK(resultRebinvaria1 == true);
}

BOOST_AUTO_TEST_SUITE_END()

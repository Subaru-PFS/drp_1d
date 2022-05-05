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
#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/mask.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <limits>
#include <vector>

#include <gsl/gsl_fit.h>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Spectrum)

BOOST_AUTO_TEST_CASE(name) {

  CSpectrum object_CSpectrum;

  const char *test_name = "toto";
  object_CSpectrum.SetName(test_name);
  BOOST_CHECK(object_CSpectrum.GetName() == test_name);
  BOOST_TEST_MESSAGE("GetName OK");
}

BOOST_AUTO_TEST_CASE(path) {

  CSpectrum object_CSpectrum;

  const char *test_path = "chemin";
  object_CSpectrum.SetFullPath(test_path);
  BOOST_CHECK(object_CSpectrum.GetFullPath() == test_path);
  BOOST_TEST_MESSAGE("GetFullPath OK");
}

BOOST_AUTO_TEST_CASE(invert) {

  CSpectrum object_CSpectrum;

  BOOST_CHECK(object_CSpectrum.InvertFlux() == true);
  BOOST_TEST_MESSAGE("InvertFlux OK");
}

BOOST_AUTO_TEST_CASE(Calcul) {

  //--------------------//
  // constructor

  CSpectrum object_CSpectrum;
  BOOST_TEST_MESSAGE("index:" << object_CSpectrum.GetSampleCount());
  BOOST_CHECK(object_CSpectrum.GetSampleCount() == 0);

  int nbmin = 0;
  int nbmax = 12;

  CSpectrumFluxAxis m_FluxAxis(nbmax);
  CSpectrumSpectralAxis m_SpectralAxis(nbmax, false);

  for (int i = nbmin; i < nbmax; ++i) {
    m_SpectralAxis[i] = i + 1;

    if (i < 5) {
      m_FluxAxis[i] = 0.0;
      m_FluxAxis.GetError()[i] = 0.0;
    } else if (i == 7) {
      m_FluxAxis[i] = std::nan("1");
      m_FluxAxis.GetError()[i] = std::nan("2");
    } else if (i == 9) {
      m_FluxAxis[i] = std::numeric_limits<double>::infinity();
      m_FluxAxis.GetError()[i] = std::numeric_limits<double>::infinity();
    } else {
      m_FluxAxis[i] = i + 2;
      m_FluxAxis.GetError()[i] = 1e-12;
    }

    BOOST_TEST_MESSAGE("m_SpectralAxis[i]:" << as_const(m_SpectralAxis)[i]);
  }

  BOOST_TEST_MESSAGE("index1:" << m_SpectralAxis.GetSamplesCount());
  BOOST_TEST_MESSAGE("index2:" << m_FluxAxis.GetSamplesCount());

  object_CSpectrum.SetSpectralAndFluxAxes(m_SpectralAxis, m_FluxAxis);

  BOOST_TEST_MESSAGE("index21:" << object_CSpectrum.GetSampleCount());
  BOOST_CHECK(object_CSpectrum.GetSampleCount() == nbmax);

  const CSpectrumFluxAxis &const_FluxAxis = m_FluxAxis;
  const CSpectrumSpectralAxis &const_SpectralAxis = m_SpectralAxis;
  const CSpectrumNoiseAxis &const_noiseAxis = m_FluxAxis.GetError();

  CSpectrum object_CSpectrum2;
  object_CSpectrum2 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum2.GetSampleCount() == nbmax);

  BOOST_TEST_MESSAGE("index22:" << object_CSpectrum2.GetSampleCount());
  BOOST_TEST_MESSAGE(
      "index23:" << as_const(object_CSpectrum2.GetFluxAxis())[0]);
  BOOST_TEST_MESSAGE(
      "index24:" << as_const(object_CSpectrum2.GetSpectralAxis())[0]);

  TFloat64List mask(nbmax, 1.);
  BOOST_TEST_MESSAGE("index31:" << mask.size());
  BOOST_TEST_MESSAGE("index32:" << mask[0]);
  BOOST_TEST_MESSAGE("index33:" << mask[nbmax - 1]);

  CSpectrum object_CSpectrum3(object_CSpectrum2, mask);
  const CSpectrumFluxAxis &const_FluxAxis3 = object_CSpectrum3.GetFluxAxis();
  const CSpectrumSpectralAxis &const_SpectralAxis3 =
      object_CSpectrum3.GetSpectralAxis();
  const CSpectrumNoiseAxis &const_noiseAxis3 =
      object_CSpectrum3.GetFluxAxis().GetError();
  //
  BOOST_CHECK(object_CSpectrum3.GetSampleCount() == nbmax);
  BOOST_CHECK(const_FluxAxis3[0] == const_FluxAxis[0]);
  BOOST_CHECK(const_noiseAxis3[0] == const_noiseAxis[0]);
  BOOST_CHECK(const_SpectralAxis3[0] == const_SpectralAxis[0]);
  BOOST_CHECK(const_FluxAxis3[nbmax - 1] == const_FluxAxis[nbmax - 1]);
  BOOST_CHECK(const_noiseAxis3[nbmax - 1] == const_noiseAxis[nbmax - 1]);
  BOOST_CHECK(const_SpectralAxis3[nbmax - 1] == const_SpectralAxis[nbmax - 1]);

  BOOST_TEST_MESSAGE("index42:" << object_CSpectrum3.GetSampleCount());
  BOOST_TEST_MESSAGE("index43:" << const_FluxAxis3[0]);
  BOOST_TEST_MESSAGE("index44:" << const_SpectralAxis3[0]);

  TFloat64List mask0(nbmax, 0.);
  CSpectrum object_CSpectrum3_bis(object_CSpectrum2, mask0);
  BOOST_CHECK(object_CSpectrum3_bis.GetSampleCount() == 0);
  BOOST_CHECK(object_CSpectrum3_bis.GetFluxAxis().isEmpty() == true);
  BOOST_CHECK(object_CSpectrum3_bis.GetFluxAxis().GetError().isEmpty() == true);
  BOOST_CHECK(object_CSpectrum3_bis.GetSpectralAxis().isEmpty() == true);

  TFloat64List mask_even;
  for (int i = nbmin; i < nbmax; ++i) {
    if (i % 2 == 0)
      mask_even.push_back(1.);
    else
      mask_even.push_back(0.);
  }
  CSpectrum object_CSpectrum3_ter(object_CSpectrum2, mask_even);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis3_ter =
      object_CSpectrum3_ter.GetFluxAxis();
  const CSpectrumSpectralAxis &const_SpectralAxis3_ter =
      object_CSpectrum3_ter.GetSpectralAxis();
  const CSpectrumNoiseAxis &const_noiseAxis3_ter =
      object_CSpectrum3_ter.GetFluxAxis().GetError();

  BOOST_CHECK(object_CSpectrum3_ter.GetSampleCount() == nbmax / 2);
  for (int i = nbmin; i < nbmax / 2; ++i) {
    BOOST_CHECK(const_FluxAxis3_ter[i] == const_FluxAxis[2 * i]);
    BOOST_CHECK(const_noiseAxis3_ter[i] == const_noiseAxis[2 * i]);
    BOOST_CHECK(const_SpectralAxis3_ter[i] == const_SpectralAxis[2 * i]);
  }

  BOOST_TEST_MESSAGE("test constructeur OK");

  //--------------------//
  // test GetMeanAndStdFluxInRange

  TFloat64Range object_CRange(1., 5.);
  TFloat64Range object_CRangeb(-1., 5.);
  TFloat64Range object_CRangec(1., 15.);

  BOOST_TEST_MESSAGE("index51:" << object_CRange.GetBegin());
  BOOST_TEST_MESSAGE("index52:" << object_CRange.GetEnd());

  Float64 mean = 10.0;
  Float64 std = 10.0;

  bool result1 =
      object_CSpectrum2.GetMeanAndStdFluxInRange(object_CRange, mean, std);
  bool result1b =
      object_CSpectrum2.GetMeanAndStdFluxInRange(object_CRangeb, mean, std);
  bool result1c =
      object_CSpectrum2.GetMeanAndStdFluxInRange(object_CRangec, mean, std);

  BOOST_TEST_MESSAGE("index53:" << result1);
  BOOST_TEST_MESSAGE("index53b:" << result1b);
  BOOST_TEST_MESSAGE("index53c:" << result1c);
  BOOST_TEST_MESSAGE(
      "index54:"
      << object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetBegin());
  BOOST_TEST_MESSAGE(
      "index55:"
      << object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetEnd());

  BOOST_CHECK(result1 == true);
  BOOST_CHECK(result1b == false);
  BOOST_CHECK(result1c == false);

  //--------------------//
  // test GetLinearRegInRange

  bool result11 =
      object_CSpectrum2.GetLinearRegInRange((object_CRange), mean, std);
  bool result11b =
      object_CSpectrum2.GetLinearRegInRange((object_CRangeb), mean, std);
  bool result11c =
      object_CSpectrum2.GetLinearRegInRange((object_CRangec), mean, std);

  BOOST_CHECK(result11 == true);
  BOOST_CHECK(result11b == false);
  BOOST_CHECK(result11c == false);

  //--------------------//

  CSpectrumFluxAxis _FluxAxis2(nbmax);
  CSpectrumFluxAxis _FluxAxis3(nbmax);
  CSpectrumFluxAxis _FluxAxis4(nbmax);
  CSpectrumFluxAxis _FluxAxis5(nbmax);
  CSpectrumFluxAxis _FluxAxis6(nbmax);

  for (int i = nbmin; i < nbmax; ++i) {
    _FluxAxis2[i] = (i + 2) * 1e+3;
    _FluxAxis3[i] = _FluxAxis2[i];
    _FluxAxis4[i] = 0.0;
    _FluxAxis6[i] = _FluxAxis2[i];
    _FluxAxis2.GetError()[i] = 1e-5;
    _FluxAxis3.GetError()[i] = 0.0;
    _FluxAxis4.GetError()[i] = _FluxAxis2.GetError()[i];
    _FluxAxis5.GetError()[i] = _FluxAxis2.GetError()[i];

    if (i < 5) {
      _FluxAxis5[i] = std::nan("5");
      _FluxAxis6.GetError()[i] = std::nan("6");
    } else if (i == 5) {
      _FluxAxis5[i] = 1e+3;
      _FluxAxis6.GetError()[i] = 1e-9;
    } else {
      _FluxAxis5[i] = std::numeric_limits<double>::infinity();
      _FluxAxis6.GetError()[i] = std::numeric_limits<double>::infinity();
    }
  }

  CSpectrum object_CSpectrum4(m_SpectralAxis, _FluxAxis2);
  CSpectrum object_CSpectrum5(m_SpectralAxis, _FluxAxis3);
  CSpectrum object_CSpectrum6(m_SpectralAxis, _FluxAxis4);
  CSpectrum object_CSpectrum7(m_SpectralAxis, _FluxAxis5);
  CSpectrum object_CSpectrum8(m_SpectralAxis, _FluxAxis6);

  //--------------------//
  // test IsFluxValid

  BOOST_CHECK(object_CSpectrum2.IsFluxValid(1, 11.1) ==
              false); // cas dans tout l'intervalle
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(1, 5.1) ==
              false); // cas dans l'intervalle 1 à 5 avec 0.0
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(1, 6.1) ==
              true); // cas dans l'intervalle 1 à 6
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(6, 11.1) ==
              false); // cas dans l'intervalle 6 à 11 avec nan et inf
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(6, 7.1) ==
              true); // cas dans l'intervalle 6 à 7
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(7, 9.1) ==
              false); // cas dans l'intervalle 7 à 9 avec nan
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(9, 9.1) ==
              true); // cas où l'intervalle est un point
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(9, 11.1) ==
              false); // cas dans l'intervalle 9 à 11 avec inf
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(10, 10.1) ==
              false); // cas où l'intervalle est un point inf
  BOOST_CHECK(object_CSpectrum2.IsFluxValid(11, 14.1) ==
              false); // cas où l'intervalle est à l'extérieur

  //--------------------//
  // test IsNoiseValid

  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(1, 11.1) ==
              false); // cas dans tout l'intervalle
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(1, 5.1) ==
              false); // cas dans l'intervalle 1 à 5 avec 0.0
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(1, 6.1) ==
              false); // cas dans l'intervalle 1 à 6
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(6, 11.1) ==
              false); // cas dans l'intervalle 6 à 11 avec nan et inf
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(6, 7.1) ==
              true); // cas dans l'intervalle 6 à 7
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(7, 9.1) ==
              false); // cas dans l'intervalle 7 à 9 avec nan
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(9, 9.1) ==
              true); // cas où l'intervalle est un point
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(9, 11.1) ==
              false); // cas dans l'intervalle 9 à 11 avec inf
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(10, 10.1) ==
              false); // cas où l'intervalle est un point inf
  BOOST_CHECK(object_CSpectrum2.IsNoiseValid(11, 14.1) ==
              false); // cas où l'intervalle est à l'extérieur

  //--------------------//
  // test correctSpectrum

  // cas où toutes les valeurs du flux et de l'erreur sont valides
  BOOST_CHECK(object_CSpectrum4.correctSpectrum(1, 11.2) == false);
  // cas où toutes les valeurs de l'erreur sont nulles, et le flux valide
  BOOST_CHECK_THROW(object_CSpectrum5.correctSpectrum(1, 11.2),
                    GlobalException);
  // cas où toutes les valeurs du flux sont nulles, et l'erreur valide
  BOOST_CHECK(object_CSpectrum6.correctSpectrum(1, 11.2) == false);

  // cas où une seule valeur du flux est valide, et l'erreur valide
  BOOST_CHECK(object_CSpectrum7.correctSpectrum(0.7, 11.3) == true);
  const TFloat64List f7 = object_CSpectrum7.GetFluxAxis().GetSamplesVector();
  const TFloat64List f7c{
      1e+2, 1e+2, 1e+2, 1e+2, 1e+2, 1e+3,
      1e+2, 1e+2, 1e+2, 1e+2, 1e+2, std::numeric_limits<double>::infinity()};
  BOOST_CHECK(f7 == f7c);
  const TFloat64List e7 =
      object_CSpectrum7.GetFluxAxis().GetError().GetSamplesVector();
  const TFloat64List e7c{1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5,
                         1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5};
  BOOST_CHECK(e7 == e7c);

  // cas où une seule valeur de l'erreur est valide, et le flux valide
  BOOST_CHECK(object_CSpectrum8.correctSpectrum(0.8, 11.3) == true);
  const TFloat64List f8 = object_CSpectrum8.GetFluxAxis().GetSamplesVector();
  const TFloat64List f8c = {7.0 * 1e+2, 7.0 * 1e+2, 7.0 * 1e+2, 7.0 * 1e+2,
                            7.0 * 1e+2, 7.0 * 1e+3, 7.0 * 1e+2, 7.0 * 1e+2,
                            7.0 * 1e+2, 7.0 * 1e+2, 7.0 * 1e+2, 13.0 * 1e+3};
  BOOST_CHECK(f8 == f8c);
  const TFloat64List e8 =
      object_CSpectrum8.GetFluxAxis().GetError().GetSamplesVector();
  const TFloat64List e8c{
      1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-9,
      1e-8, 1e-8, 1e-8, 1e-8, 1e-8, std::numeric_limits<double>::infinity()};
  BOOST_CHECK(e8 == e8c);

  // cas dans l'intervalle 1 à 5 avec 0.0
  CSpectrum object_CSpectrum2copy = object_CSpectrum2;
  BOOST_CHECK_THROW(object_CSpectrum2copy.correctSpectrum(1, 5.4),
                    GlobalException);

  // cas dans l'intervalle 1 à 6
  CSpectrum object_CSpectrum9 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum9.correctSpectrum(1, 6.4) == true);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis9 = object_CSpectrum9.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis9 =
      object_CSpectrum9.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i < 5) { // valeurs corrigées
      BOOST_CHECK(const_FluxAxis9[i] == 0.7);
      BOOST_CHECK(const_noiseAxis9[i] == 1e-11);
    } else if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis9[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis9[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis9[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis9[i] == const_noiseAxis[i]);
    }
  }

  // cas dans l'intervalle 6 à 7
  CSpectrum object_CSpectrum10 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum10.correctSpectrum(6, 7.4) == false);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis10 = object_CSpectrum10.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis10 =
      object_CSpectrum10.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis10[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis10[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis10[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis10[i] == const_noiseAxis[i]);
    }
  }

  // cas dans l'intervalle 7 à 9 avec nan
  CSpectrum object_CSpectrum11 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum11.correctSpectrum(7, 9.4) == true);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis11 = object_CSpectrum11.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis11 =
      object_CSpectrum11.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 7) { // nan en position 8 : valeur corrigée
      BOOST_CHECK(const_FluxAxis11[i] == 0.8);
      BOOST_CHECK(const_noiseAxis11[i] == 1e-11);
    } else {
      BOOST_CHECK(const_FluxAxis11[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis11[i] == const_noiseAxis[i]);
    }
  }

  // cas où l'intervalle est un point
  CSpectrum object_CSpectrum12 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum12.correctSpectrum(9, 9.4) == false);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis12 = object_CSpectrum12.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis12 =
      object_CSpectrum12.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis12[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis12[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis12[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis12[i] == const_noiseAxis[i]);
    }
  }

  // cas dans l'intervalle 9 à 11 avec inf
  CSpectrum object_CSpectrum13 = object_CSpectrum;
  BOOST_CHECK(object_CSpectrum13.correctSpectrum(9, 11.4) == true);
  // const access
  const CSpectrumFluxAxis &const_FluxAxis13 = object_CSpectrum13.GetFluxAxis();
  const CSpectrumNoiseAxis &const_noiseAxis13 =
      object_CSpectrum13.GetFluxAxis().GetError();

  for (int i = nbmin; i < nbmax; ++i) {
    if (i == 9) { // inf en position 10 : valeur corrigée
      BOOST_CHECK(const_FluxAxis13[i] == 1.0);
      BOOST_CHECK(const_noiseAxis13[i] == 1e-11);
    } else if (i == 7) { // nan en position 8
      BOOST_CHECK(std::isnan(const_FluxAxis13[i]) == true);
      BOOST_CHECK(std::isnan(const_noiseAxis13[i]) == true);
    } else {
      BOOST_CHECK(const_FluxAxis13[i] == const_FluxAxis[i]);
      BOOST_CHECK(const_noiseAxis13[i] == const_noiseAxis[i]);
    }
  }

  // cas après correction où l'intervalle est un point inf
  BOOST_CHECK(object_CSpectrum13.correctSpectrum(10, 10.5) == false);

  // cas après correction dans l'intervalle 6 à 9 avec nan
  BOOST_CHECK(object_CSpectrum11.correctSpectrum(6, 9.5) == false);

  // cas après correction dans l'intervalle 1 à 7
  BOOST_CHECK(object_CSpectrum9.correctSpectrum(1, 7.5) == false);

  // cas où l'intervalle est à l'extérieur
  BOOST_CHECK(object_CSpectrum2.correctSpectrum(11, 14.5) == false);

  //--------------------//
  // test GetLambdaRange

  Float64 intervalle =
      object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetLength();

  BOOST_CHECK_CLOSE(intervalle, object_CSpectrum2.GetLambdaRange().GetLength(),
                    1e-12);

  BOOST_TEST_MESSAGE(
      "result9A:"
      << object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetLength());
  BOOST_TEST_MESSAGE(
      "index55A:"
      << object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetBegin());
  BOOST_TEST_MESSAGE(
      "index55B:"
      << object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetEnd());

  //--------------------//
  // test ConvertToLinearScale

  bool result2 = object_CSpectrum2.ConvertToLinearScale();
  CSpectrumSpectralAxis wavAxis1 = object_CSpectrum2.GetSpectralAxis();
  bool result3 = wavAxis1.ConvertToLinearScale();

  BOOST_CHECK(result2 == result3);

  BOOST_TEST_MESSAGE("result2:" << result2);
  BOOST_TEST_MESSAGE("result3:" << result3);

  //--------------------//
  /// test ConvertToLogScale

  bool result4 = object_CSpectrum2.ConvertToLogScale();
  CSpectrumSpectralAxis wavAxis2 = object_CSpectrum2.GetSpectralAxis();
  bool result5 = wavAxis2.ConvertToLogScale();

  BOOST_CHECK(result4 == result5);

  BOOST_TEST_MESSAGE("result4:" << result4);
  BOOST_TEST_MESSAGE("result5:" << result5);

  //--------------------//
  // test GetMeanResolution

  Float64 result6 = object_CSpectrum2.GetMeanResolution();
  Float64 result7 = object_CSpectrum2.GetSpectralAxis().GetMeanResolution();

  BOOST_CHECK_CLOSE(result6, result7, 1e-12);

  BOOST_TEST_MESSAGE("result6:" << result6);
  BOOST_TEST_MESSAGE("result7:" << result7);

  //--------------------//
  // test GetResolution

  Float64 result8 = object_CSpectrum2.GetResolution();
  Float64 result9 = object_CSpectrum2.GetSpectralAxis().GetResolution();

  BOOST_CHECK_CLOSE(result8, result9, 1e-12);

  BOOST_TEST_MESSAGE("result8:" << result8);
  BOOST_TEST_MESSAGE("result9:" << result9);

  //--------------------//
  // //test removeContinuum

  CContinuumIrregularSamplingMedian remover2;

  for (int i = nbmin; i < nbmax; ++i) {
    m_FluxAxis[i] = 2.0;
  }
  object_CSpectrum2.SetSpectralAndFluxAxes(m_SpectralAxis, m_FluxAxis);

  BOOST_CHECK(object_CSpectrum2.RemoveContinuum(remover2) == true);
  BOOST_TEST_MESSAGE(
      "test Remove:" << object_CSpectrum2.RemoveContinuum(remover2));
}

BOOST_AUTO_TEST_CASE(ExtractTest) {
  const CSpectrumAxis axis({1., 2., 3., 4., 5.});
  const CSpectrumSpectralAxis spcAxis({1., 2., 3., 4., 5.});
  const CSpectrumNoiseAxis noiseAxis({-1., -2., -3., -4., -5.});
  const CSpectrumFluxAxis fluxAxis(CSpectrumAxis({2., 4., 6., 8., 10.}),
                                   noiseAxis);
  const CSpectrum spc(spcAxis, fluxAxis);

  Int32 istart = 1;
  Int32 iend = 3;
  Int32 s = iend - istart + 1;

  const TFloat64List correctSpcAxis{2., 3., 4.};
  const TFloat64List correctFluxAxis{4., 6., 8};
  const TFloat64List correctNoiseAxis{-2., -3., -4.};

  ///////////////////////////
  const CSpectrumAxis axis_extract = axis.extract(istart, iend);
  const CSpectrumSpectralAxis spcAxis_extract = spcAxis.extract(istart, iend);
  const CSpectrumNoiseAxis noiseAxis_extract = noiseAxis.extract(istart, iend);
  const CSpectrumFluxAxis fluxAxis_extract = fluxAxis.extract(istart, iend);

  const TFloat64List &extractedAxis1 = axis_extract.GetSamplesVector();
  const TFloat64List &extractedSpcAxis1 = spcAxis_extract.GetSamplesVector();
  const TFloat64List &extractedFluxAxis1 = fluxAxis_extract.GetSamplesVector();
  const TFloat64List &extractedNoiseAxis1 =
      noiseAxis_extract.GetSamplesVector();

  BOOST_CHECK(extractedAxis1.size() == s);
  BOOST_CHECK(extractedSpcAxis1.size() == s);
  BOOST_CHECK(extractedFluxAxis1.size() == s);
  BOOST_CHECK(extractedNoiseAxis1.size() == s);
  ////////////////////////////
  // read results
  const CSpectrum extractedSpc = spc.extract(istart, iend);
  const TFloat64List &extractedSpcAxis =
      extractedSpc.GetSpectralAxis().GetSamplesVector();
  const TFloat64List &extractedFluxAxis =
      extractedSpc.GetFluxAxis().GetSamplesVector();
  const TFloat64List &extractedNoiseAxis =
      extractedSpc.GetFluxAxis().GetError().GetSamplesVector();

  // check results
  // check size is correct
  BOOST_CHECK(extractedSpcAxis.size() == s);
  BOOST_CHECK(extractedFluxAxis.size() == s);
  BOOST_CHECK(extractedNoiseAxis.size() == s);

  BOOST_CHECK_EQUAL_COLLECTIONS(extractedSpcAxis.begin(),
                                extractedSpcAxis.end(), correctSpcAxis.begin(),
                                correctSpcAxis.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(extractedFluxAxis.begin(),
                                extractedFluxAxis.end(),
                                correctFluxAxis.begin(), correctFluxAxis.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(
      extractedNoiseAxis.begin(), extractedNoiseAxis.end(),
      correctNoiseAxis.begin(), correctNoiseAxis.end());
}
BOOST_AUTO_TEST_SUITE_END()

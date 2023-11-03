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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "tests/src/tool/inputContextLight.h"

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace NSEpic;
using namespace std;

namespace bfs = boost::filesystem;

// //-----------------------------------------------------------------------------

const std::string jsonString =
    "{\"lsf\" : {\"lsfType\" : \"gaussianConstantWidth\" , \"width\" : 0.04},"
    "\"galaxy\" : {\"method\" : \"templateFittingSolver\", "
    "\"templateFittingSolver\" : {\"dustfit\" : true, \"extinction\" : true}}, "
    "\"qso\" : {\"method\" : \"templateFittingSolver\", "
    "\"templateFittingSolver\" : {\"dustfit\" : true, \"extinction\" : true}}, "
    "\"templateCatalog\" : {\"continuumRemoval\" : "
    "{ \"medianKernelWidth\" : 40.0, "
    "\"medianEvenReflection\" : false, "
    "\"method\" : \"irregularSamplingMedian\"}}}";

class fixture_MeiskinCorrectionTest {
public:
  TScopeStack scopeStack;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  std::shared_ptr<CLSF> LSFConstantWidth =
      fixture_LSFGaussianConstantWidth(scopeStack, jsonString).LSF;
  std::shared_ptr<CLSF> LSFConstantResol =
      fixture_LSFGaussianConstantResolution(scopeStack).LSF;
  TFloat64List z_bins = fixture_MeiskinCorrection().z_bins;
  Int32 idxCount = fixture_MeiskinCorrection().idxCount;
};

BOOST_FIXTURE_TEST_SUITE(fluxcorrectionmeiksin_test,
                         fixture_MeiskinCorrectionTest)

BOOST_AUTO_TEST_CASE(overall_test) {
  Int32 z_index = igmCorrectionMeiksin->getRedshiftIndex(0.5);
  BOOST_CHECK(z_index == -1);
  z_index = igmCorrectionMeiksin->getRedshiftIndex(2.2);
  BOOST_CHECK(z_index == 0);
  z_index = igmCorrectionMeiksin->getRedshiftIndex(2.7);
  BOOST_CHECK(z_index == 1);
  z_index = igmCorrectionMeiksin->getRedshiftIndex(3.0);
  BOOST_CHECK(z_index == 1);

  // test .h methods
  BOOST_CHECK(igmCorrectionMeiksin->isConvolved() == false);
  BOOST_CHECK(igmCorrectionMeiksin->getIdxCount() == idxCount);
  BOOST_CHECK(igmCorrectionMeiksin->getLambdaMin() ==
              igmCorrectionMeiksin->m_rawCorrections[0].lbda.front());
  BOOST_CHECK(igmCorrectionMeiksin->getLambdaMax() == RESTLAMBDA_LYA);
  BOOST_CHECK(igmCorrectionMeiksin->getRedshiftBins() == z_bins);

  // rawGrid & fineGrid has same size for now
  Float64 finegridstep = 0.05;
  TFloat64Range range(igmCorrectionMeiksin->m_rawCorrections[0].lbda.front(),
                      igmCorrectionMeiksin->m_rawCorrections[0].lbda.back());
  TFloat64List finelbdaGrid = range.SpreadOver(finegridstep);

  TFloat64Range zbin(z_bins[1], z_bins[2]);
  igmCorrectionMeiksin->m_convolRange = TFloat64Range(3900, 12500);
  igmCorrectionMeiksin->m_fineLambdaSize = finelbdaGrid.size();
  //
  TFloat64List convolvedArray = igmCorrectionMeiksin->ConvolveByLSFOneCurve(
      igmCorrectionMeiksin->m_rawCorrections[0].fluxcorr[0],
      igmCorrectionMeiksin->m_rawCorrections[0].lbda, finelbdaGrid, zbin,
      LSFConstantWidth);

  for (std::size_t i = 0; i < finelbdaGrid.size(); i++) {
    BOOST_CHECK_CLOSE(convolvedArray[i],
                      igmCorrectionMeiksin->m_rawCorrections[0].fluxcorr[0][i],
                      1e-8);
  }

  // empty array
  BOOST_CHECK_THROW(igmCorrectionMeiksin->ConvolveByLSFOneCurve(
                        TFloat64List{},
                        igmCorrectionMeiksin->m_rawCorrections[0].lbda,
                        finelbdaGrid, zbin, LSFConstantWidth),
                    GlobalException);

  // no overlapping range
  zbin.Set(0.5, 1.0);
  convolvedArray = igmCorrectionMeiksin->ConvolveByLSFOneCurve(
      igmCorrectionMeiksin->m_rawCorrections[0].fluxcorr[0],
      igmCorrectionMeiksin->m_rawCorrections[0].lbda, finelbdaGrid, zbin,
      LSFConstantWidth);
  BOOST_CHECK(
      convolvedArray ==
      TFloat64List(igmCorrectionMeiksin->m_rawCorrections[0].lbda.size(), 0.));

  // test getWaveIndex :
  Int32 w_index = igmCorrectionMeiksin->getWaveIndex(
      igmCorrectionMeiksin->m_rawCorrections[0].lbda.front());
  BOOST_CHECK(w_index == 0);

  w_index = igmCorrectionMeiksin->getWaveIndex(1048.06999999998);
  Int32 w_ref = std::round((1048.06999999998 - 1000.02999999997) / 0.05);
  BOOST_CHECK(w_index == w_ref);

  w_index = igmCorrectionMeiksin->getWaveIndex(1220.);
  BOOST_CHECK(w_index ==
              igmCorrectionMeiksin->m_rawCorrections[0].lbda.size() - 1);

  // test getWaveVector :
  TFloat64Range w_range(1213., 1213.2);
  TFloat64List w_vec = igmCorrectionMeiksin->getWaveVector(w_range);

  TFloat64List w_vec_ref = {1213.02999999983, 1213.07999999983,
                            1213.12999999983, 1213.17999999983};
  for (std::size_t i = 0; i < w_vec.size(); i++) {
    BOOST_CHECK_CLOSE(w_vec[i], w_vec_ref[i], 1e-8);
  }

  TFloat64Range w_range_raw(1213., 1213.2);
  w_vec = igmCorrectionMeiksin->getWaveVector(w_range_raw, true);
  for (std::size_t i = 0; i < w_vec.size(); i++) {
    BOOST_CHECK_CLOSE(w_vec[i], w_vec_ref[i], 1e-8);
  }

  TFloat64Range w_range_2(1000.93, 1000.94);
  w_vec = igmCorrectionMeiksin->getWaveVector(w_range_2, true);
  BOOST_CHECK(w_vec.size() == 0);
}

BOOST_AUTO_TEST_CASE(convolveByLSF_test) {
  // to test kernel 1 use width = 0.045 -> give the same thing
  // rawGrid & fineGrid has same size for now
  Float64 finegridstep = 0.05;
  TFloat64Range range(igmCorrectionMeiksin->m_rawCorrections[0].lbda.front(),
                      igmCorrectionMeiksin->m_rawCorrections[0].lbda.back());
  TFloat64List finelbdaGrid = range.SpreadOver(finegridstep);

  // test convolve with kernel = 1
  TFloat64Range lbdaRange(1000, 12500);
  TScopeStack scopeStack;
  igmCorrectionMeiksin->convolveByLSF(LSFConstantWidth, lbdaRange);

  //  range [2.0, 2.5]
  Float64 redshift = 2.1;
  // 7 flux per curve
  for (std::size_t i = 0; i < igmCorrectionMeiksin->getIdxCount(); i++) {
    for (std::size_t j = 0; j < finelbdaGrid.size(); j++) {
      Float64 corr =
          igmCorrectionMeiksin->getCorrection(redshift, i, finelbdaGrid[j]);
      if (finelbdaGrid[j] > RESTLAMBDA_LYA)
        BOOST_CHECK(corr == 1.);
      else
        BOOST_CHECK_CLOSE(
            corr, igmCorrectionMeiksin->m_rawCorrections[0].fluxcorr[0][j],
            1e-8);
    }
  }
  //  range [2.5, 3.0]
  redshift = 2.7;
  // 7 flux per curve
  for (std::size_t i = 0; i < igmCorrectionMeiksin->getIdxCount(); i++) {
    for (std::size_t j = 0; j < finelbdaGrid.size(); j++) {
      Float64 corr =
          igmCorrectionMeiksin->getCorrection(redshift, i, finelbdaGrid[j]);
      if (finelbdaGrid[j] > RESTLAMBDA_LYA)
        BOOST_CHECK(corr == 1.);
      else
        BOOST_CHECK_CLOSE(
            corr, igmCorrectionMeiksin->m_rawCorrections[1].fluxcorr[0][j],
            1e-8);
    }
  }

  // test getCorrection
  Float64 lambdaRest = 1000.8;
  Float64 corr = igmCorrectionMeiksin->getCorrection(redshift, 0, lambdaRest);
  BOOST_CHECK(corr = 0.752732320000021);

  lambdaRest = 1217.;
  corr = igmCorrectionMeiksin->getCorrection(redshift, 0, lambdaRest);
  BOOST_CHECK(corr = 1);

  Int32 meiksinIdx = -1;
  lambdaRest = 1214.;
  BOOST_CHECK_THROW(
      igmCorrectionMeiksin->getCorrection(redshift, meiksinIdx, lambdaRest),
      GlobalException);

  redshift = 1.5;
  meiksinIdx = 0;
  corr = igmCorrectionMeiksin->getCorrection(redshift, meiksinIdx, lambdaRest);
  BOOST_CHECK(corr = 1);

  // test convolve with real convolution
  redshift = 2.1;
  igmCorrectionMeiksin->convolveByLSF(LSFConstantResol, lbdaRange);
  for (std::size_t i = 0; i < finelbdaGrid.size(); i++) {
    Float64 corr = igmCorrectionMeiksin->getCorrection(redshift, meiksinIdx,
                                                       finelbdaGrid[i]);
    if (finelbdaGrid[i] > RESTLAMBDA_LYA)
      BOOST_CHECK(corr == 1.);
    else
      BOOST_CHECK_CLOSE(
          corr, igmCorrectionMeiksin->m_rawCorrections[0].fluxcorr[0][i], 1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()

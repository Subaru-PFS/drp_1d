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

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(fluxcorrectionmeiksin_test)

namespace bfs = boost::filesystem;

// //-----------------------------------------------------------------------------

const std::string jsonString =
    "{\"galaxy\" : {\"method\" : \"TemplateFittingSolve\", "
    "\"TemplateFittingSolve\" : {\"dustfit\" : true, \"extinction\" : true}}, "
    "\"qso\" : {\"method\" : \"TemplateFittingSolve\", "
    "\"TemplateFittingSolve\" : {\"dustfit\" : true, \"extinction\" : true}}, "
    "\"templateCatalog\" : {\"continuumRemoval\" : "
    "{ \"medianKernelWidth\" : 40.0, "
    "\"medianEvenReflection\" : false, "
    "\"method\" : \"IrregularSamplingMedian\"}}}";

TFloat64List lbda = {
    1213.01999999983, 1213.06999999983, 1213.11999999983, 1213.16999999983,
    1213.21999999983, 1213.26999999983, 1213.31999999983, 1213.36999999983,
    1213.41999999983, 1213.46999999983, 1213.51999999983, 1213.56999999983,
    1213.61999999983, 1213.66999999983, 1213.71999999983, 1213.76999999983,
    1213.81999999983, 1213.86999999983, 1213.91999999983, 1213.96999999983,
    1214.01999999983, 1214.06999999983, 1214.11999999983, 1214.16999999983,
    1214.21999999983, 1214.26999999983, 1214.31999999983, 1214.36999999983,
    1214.41999999983, 1214.46999999983, 1214.51999999983, 1214.56999999983,
    1214.61999999983, 1214.66999999983, 1214.71999999983, 1214.76999999983,
    1214.81999999983, 1214.86999999983, 1214.91999999983, 1214.96999999983,
    1215.01999999983, 1215.06999999983, 1215.11999999983, 1215.16999999983,
    1215.21999999983, 1215.26999999983, 1215.31999999983, 1215.36999999983,
    1215.41999999983, 1215.46999999983, 1215.51999999983, 1215.56999999983,
    1215.61999999983, 1215.66999999983, 1215.72,          1215.77,
    1215.82,          1215.87,          1215.92,          1215.97,
    1216.02,          1216.07,          1216.12,          1216.17,
    1216.22,          1216.27,          1216.32,          1216.37,
    1216.42,          1216.47,          1216.52,          1216.57,
    1216.62,          1216.67,          1216.72,          1216.77,
    1216.82,          1216.87,          1216.92,          1216.97,
    1217.02,          1217.07,          1217.12,          1217.17,
    1217.22,          1217.27,          1217.32,          1217.37,
    1217.42,          1217.47,          1217.52,          1217.57,
    1217.62,          1217.67,          1217.72,          1217.77,
    1217.82,          1217.87,          1217.92,          1217.97,
    1218.02,          1218.07,          1218.12,          1218.17,
    1218.22,          1218.27,          1218.32,          1218.37,
    1218.42,          1218.47,          1218.52,          1218.57,
    1218.62,          1218.67,          1218.72,          1218.77,
    1218.82,          1218.87,          1218.92,          1218.97,
    1219.02};
TFloat64List flux = {0.70983580837752,
                     0.70977125332827,
                     0.70970669827902,
                     0.70964214322977,
                     0.70957758818052,
                     0.70951303313127,
                     0.70944847808202,
                     0.70938392303277,
                     0.70931936798352,
                     0.70925481293427,
                     0.709190257885021,
                     0.709125702835771,
                     0.709061147786521,
                     0.708996592737271,
                     0.708932037688021,
                     0.708867482638771,
                     0.708802927589521,
                     0.708738372540271,
                     0.708673817491021,
                     0.708609262441771,
                     0.707699961319093,
                     0.705523541068195,
                     0.703347120817297,
                     0.701170700566399,
                     0.698994280315501,
                     0.696817860064603,
                     0.694641439813705,
                     0.692465019562807,
                     0.690288599311909,
                     0.688112179061011,
                     0.685935758810113,
                     0.683759338559215,
                     0.681582918308317,
                     0.679406498057419,
                     0.677230077806521,
                     0.675053657555623,
                     0.672877237304725,
                     0.670700817053827,
                     0.668524396802929,
                     0.666347976552031,
                     0.664171556301133,
                     0.661995136050235,
                     0.659818715799337,
                     0.657642295548439,
                     0.655465875297541,
                     0.653289455046643,
                     0.651113034795745,
                     0.648936614544847,
                     0.646760194293949,
                     0.644583774043051,
                     0.642407353792153,
                     0.640230933541255,
                     0.638054513290357,
                     0.635878093039459,
                     0.906211794830329,
                     0.922959688610614,
                     0.939707582390899,
                     0.956455476171183,
                     0.973203369951468,
                     0.989951263731753,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1,
                     1};

class MyInputContext {
public:
  std::shared_ptr<CParameterStore> paramStore;
  std::shared_ptr<CLSF> constantWidthLSF;
  std::shared_ptr<CLSF> constantResolutionLSF;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin;

  void InitContext() {
    TScopeStack scopeStack;
    paramStore = std::make_shared<CParameterStore>(scopeStack);
    paramStore->FromString(jsonString);

    paramStore->Set("LSF.width", 0.05);
    std::shared_ptr<TLSFArguments> args =
        std::make_shared<TLSFGaussianConstantWidthArgs>(paramStore);
    constantWidthLSF = LSFFactory.Create("GaussianConstantWidth", args);

    paramStore->Set("LSF.resolution", 4300.);
    args = std::make_shared<TLSFGaussianConstantResolutionArgs>(paramStore);
    constantResolutionLSF =
        LSFFactory.Create("GaussianConstantResolution", args);
  }

  std::shared_ptr<CParameterStore> GetParameterStore() { return paramStore; }
  std::shared_ptr<CLSF> GetConstantWidthLSF() { return constantWidthLSF; }
  std::shared_ptr<CLSF> GetConstantResolutionLSF() {
    return constantResolutionLSF;
  }
};

BOOST_AUTO_TEST_CASE(overall_test) {
  // test getRedshiftIndex
  std::vector<MeiksinCorrection> meiskinCorr;
  meiskinCorr.push_back(MeiksinCorrection(lbda, {flux}));
  TFloat64List z_bins = {3.0, 3.5};
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> correction =
      std::make_shared<CSpectrumFluxCorrectionMeiksin>(meiskinCorr, z_bins);

  Int32 z_index = correction->getRedshiftIndex(0.5);
  BOOST_CHECK(z_index == -1);
  z_index = correction->getRedshiftIndex(3.1);
  BOOST_CHECK(z_index == 0);
  z_index = correction->getRedshiftIndex(3.4);
  BOOST_CHECK(z_index == 0);
  z_index = correction->getRedshiftIndex(3.5);
  BOOST_CHECK(z_index == 0);

  // test .h methods
  BOOST_CHECK(correction->isConvolved() == false);
  BOOST_CHECK(correction->getIdxCount() == 7);
  BOOST_CHECK(correction->getLambdaMin() == lbda.front());
  BOOST_CHECK(correction->getLambdaMax() == RESTLAMBDA_LYA);
  BOOST_CHECK(correction->getRedshiftBins() == z_bins);

  // rawGrid & fineGrid has same size for now
  Float64 finegridstep = 0.05;
  TFloat64Range range(lbda[0], lbda.back());
  TFloat64List finelbdaGrid = range.SpreadOver(finegridstep);

  MyInputContext ctx;
  ctx.InitContext();

  // test ConvolveByLSFOneCurve
  TFloat64Range zbin(3.0, 3.5);
  correction->m_convolRange = TFloat64Range(3900, 12500);
  correction->m_fineLambdaSize = finelbdaGrid.size();
  TFloat64List convolvedArray = correction->ConvolveByLSFOneCurve(
      flux, lbda, finelbdaGrid, zbin, ctx.GetConstantWidthLSF());

  for (std::size_t i = 0; i < finelbdaGrid.size(); i++) {
    BOOST_CHECK(convolvedArray[i] == flux[i]);
  }

  // empty array
  BOOST_CHECK_THROW(
      correction->ConvolveByLSFOneCurve(TFloat64List{}, lbda, finelbdaGrid,
                                        zbin, ctx.GetConstantWidthLSF()),
      GlobalException);

  // no overlapping range
  zbin.Set(0.5, 1.0);
  convolvedArray = correction->ConvolveByLSFOneCurve(
      flux, lbda, finelbdaGrid, zbin, ctx.GetConstantWidthLSF());
  BOOST_CHECK(convolvedArray == TFloat64List(lbda.size(), 0.));

  // test getWaveIndex :
  Int32 w_index = correction->getWaveIndex(1213.01999999983);
  BOOST_CHECK(w_index == 0);

  w_index = correction->getWaveIndex(1215.1);
  Int32 w_ref = std::round((1215.1 - 1213.01999999983) / 0.05);
  BOOST_CHECK(w_index == w_ref);

  w_index = correction->getWaveIndex(1220.);
  BOOST_CHECK(w_index == lbda.size() - 1);

  // test getWaveVector :
  TFloat64Range w_range(1213., 1213.2);
  TFloat64List w_vec = correction->getWaveVector(w_range);

  TFloat64List w_vec_ref = {1213.01999999983, 1213.06999999983,
                            1213.11999999983, 1213.16999999983};
  for (std::size_t i = 0; i < w_vec.size(); i++) {
    BOOST_CHECK_CLOSE(w_vec[i], w_vec_ref[i], 1e-8);
  }

  TFloat64Range w_range_raw(1213., 1213.2);
  w_vec = correction->getWaveVector(w_range_raw, true);
  for (std::size_t i = 0; i < w_vec.size(); i++) {
    BOOST_CHECK_CLOSE(w_vec[i], w_vec_ref[i], 1e-8);
  }

  TFloat64Range w_range_2(1213.14, 1213.15);
  w_vec = correction->getWaveVector(w_range_2, true);
  BOOST_CHECK(w_vec.size() == 0);
}

BOOST_AUTO_TEST_CASE(convolveByLSF_test) {
  MyInputContext ctx;
  ctx.InitContext();

  std::vector<MeiksinCorrection> meiskinCorr;
  meiskinCorr.push_back(MeiksinCorrection(lbda, {flux}));
  TFloat64List z_bins = {3.0, 3.5};
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> correction =
      std::make_shared<CSpectrumFluxCorrectionMeiksin>(meiskinCorr, z_bins);

  // to test kernel 1 use width = 0.05 -> give the same thing
  // rawGrid & fineGrid has same size for now
  Float64 finegridstep = 0.05;
  TFloat64Range range(lbda[0], lbda.back());
  TFloat64List finelbdaGrid = range.SpreadOver(finegridstep);

  // test convolve with kernel = 1
  Float64 redshift = 3.5;
  Int32 meiksinIdx = 0;
  Float64 lambdaRest = 1213.1;
  TFloat64Range lbdaRange(3900, 12500);
  correction->convolveByLSF(ctx.GetConstantWidthLSF(), lbdaRange);

  for (std::size_t i = 0; i < finelbdaGrid.size(); i++) {
    Float64 corr =
        correction->getCorrection(redshift, meiksinIdx, finelbdaGrid[i]);
    if (finelbdaGrid[i] > RESTLAMBDA_LYA)
      BOOST_CHECK(corr == 1.);
    else
      BOOST_CHECK(corr == flux[i]);
  }

  // test getCorrection
  Float64 corr = correction->getCorrection(redshift, meiksinIdx, lambdaRest);
  BOOST_CHECK(corr = 0.70977125332827);

  lambdaRest = 1217.;
  corr = correction->getCorrection(redshift, meiksinIdx, lambdaRest);
  BOOST_CHECK(corr = 1);

  meiksinIdx = -1;
  lambdaRest = 1214.;
  BOOST_CHECK_THROW(correction->getCorrection(redshift, meiksinIdx, lambdaRest),
                    GlobalException);

  redshift = 1.5;
  meiksinIdx = 0;
  corr = correction->getCorrection(redshift, meiksinIdx, lambdaRest);
  BOOST_CHECK(corr = 1);

  // test convolve with real convolution
  correction->convolveByLSF(ctx.GetConstantResolutionLSF(), lbdaRange);
  for (std::size_t i = 0; i < finelbdaGrid.size(); i++) {
    Float64 corr =
        correction->getCorrection(redshift, meiksinIdx, finelbdaGrid[i]);
    if (finelbdaGrid[i] > RESTLAMBDA_LYA)
      BOOST_CHECK(corr == 1.);
    else
      BOOST_CHECK(corr != flux[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()

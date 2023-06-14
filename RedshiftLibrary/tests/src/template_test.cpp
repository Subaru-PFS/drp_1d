
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
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "tests/src/tool/inputContextLight.h"

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <cmath>

using namespace NSEpic;

class fixture_TemplateTest {
public:
  fixture_TemplateTest() {
    igmCorrectionMeiksin->convolveByLSF(LSF,
                                        fixture_MeiskinCorrection().lbdaRange);
  }
  TScopeStack scopeStack;
  // CTemplate tplStar = fixture_TemplateStarLight().tplStarLight;
  CTemplate tplStar = fixture_TemplateStar().tplStar;
  // std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
  //     fixture_CalzettiCorrectionLight().ismCorrectionCalzetti;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  // std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
  //     fixture_MeiskinCorrectionLight().igmCorrectionMeiksin;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  CSpectrumSpectralAxis spcAxis =
      fixture_TemplateStar().tplStar.GetSpectralAxis();
  CSpectrumFluxAxis fluxAxis = fixture_TemplateStar().tplStar.GetRawFluxAxis();
  std::shared_ptr<CLSF> LSF = fixture_LSFGaussianConstantWidth(scopeStack).LSF;

  // list
  TFloat64List spcAxisList = fixture_TemplateStar().spcAxisList;
  TFloat64List fluxAxisList = fixture_TemplateStar().fluxAxisList;

  // size
  Int32 spcAxisSize = fixture_TemplateStar().spcAxisSize;
  Int32 idxCount = fixture_MeiskinCorrection().idxCount;
};

BOOST_FIXTURE_TEST_SUITE(Template, fixture_TemplateTest)

BOOST_AUTO_TEST_CASE(Constructor_test) {
  // create template
  CTemplate tpl;
  BOOST_CHECK(tpl.GetCategory() == "");
  BOOST_CHECK(tpl.m_Name == "");
  BOOST_CHECK(tpl.GetSampleCount() == 0);

  CTemplate tpl2("name", "category");
  BOOST_CHECK(tpl2.m_Name == "name");
  BOOST_CHECK(tpl2.GetCategory() == "category");
  BOOST_CHECK(tpl2.GetSampleCount() == 0);

  BOOST_CHECK(tplStar.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  // copy and copy assignement
  CTemplate tpl4(tplStar);
  BOOST_CHECK(tpl4.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tpl4.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  CTemplate tpl5(tplStar);
  TAxisSampleList rawFlux = tpl5.GetRawFluxAxis_().GetSamplesVector();
  rawFlux.pop_back();
  tpl5.GetRawFluxAxis_().setSamplesVector(rawFlux);
  BOOST_CHECK_THROW(CTemplate tpl6(tpl5), GlobalException);

  CTemplate tpl7;
  tpl7 = tplStar;
  BOOST_CHECK(tpl7.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tpl7.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  // move and move assignement
  CTemplate tpl8(std::move(tplStar));
  BOOST_CHECK(tpl8.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tpl8.GetFluxAxis().GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tplStar.GetSampleCount() == 0);

  BOOST_CHECK_THROW(CTemplate tpl8b(std::move(tpl5)), GlobalException);

  CTemplate tpl9;
  tpl9 = std::move(tpl8);
  BOOST_CHECK(tpl9.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tpl9.GetFluxAxis().GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tpl8.GetSampleCount() == 0);

  // InitIsmIgmConfig
  TFloat64List maskList(spcAxisSize, 1);
  maskList[0] = 0;
  tpl9.InitIsmIgmConfig(0, 40, 2.8, ismCorrectionCalzetti,
                        igmCorrectionMeiksin);
  CTemplate tpl10(tpl9, maskList);
  BOOST_CHECK(tpl10.m_NoIsmIgmFluxAxis.GetSamplesCount() == spcAxisSize - 1);
  BOOST_CHECK(tpl10.m_NoIsmIgmFluxAxis[0] == fluxAxisList[1]);
  BOOST_CHECK(tpl10.m_computedDustCoeff.size() == spcAxisSize - 1);
  BOOST_CHECK(tpl10.m_computedDustCoeff[0] == 1.);
  BOOST_CHECK(tpl10.m_computedMeiksingCoeff.size() == spcAxisSize - 1);
  BOOST_CHECK(tpl10.m_computedMeiksingCoeff[0] == 1.);

  TFloat64List maskList2(spcAxisSize, 0);
  tpl7.InitIsmIgmConfig(0, 2, 2.8, ismCorrectionCalzetti, igmCorrectionMeiksin);
  maskList2[3] = 1.;
  CTemplate tpl11(tpl7, maskList2);
  BOOST_CHECK(tpl11.m_NoIsmIgmFluxAxis.GetSamplesCount() == 0);
  BOOST_CHECK(tpl11.m_computedDustCoeff.size() == 0);
  BOOST_CHECK(tpl11.m_computedMeiksingCoeff.size() == 0);
}

BOOST_AUTO_TEST_CASE(Save) {
  TFloat64List array = {0., 2., 3., 6.};
  TFloat64List flux = {0.1, 0.2, 0.3, 0.4};
  CSpectrumSpectralAxis spectralAxis(array);
  CSpectrumFluxAxis fluxAxis(flux);
  CTemplate tmpl("name", "category", spectralAxis, fluxAxis);

  // linear
  boost::filesystem::path tempfile =
      boost::filesystem::unique_path("tst_%%%%%%%%%%.txt");
  const char *filename = tempfile.c_str();
  BOOST_CHECK_NO_THROW(tmpl.Save(filename));
  boost::filesystem::remove(filename);

  // bad filename
  BOOST_CHECK(tmpl.Save("") == false);
}

BOOST_AUTO_TEST_CASE(InitIsmIgmConfig_test) {

  // InitIsmIgmConfig with kstart, kend & redshift
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86),
                    GlobalException);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(spcAxisSize, spcAxisSize + 2, 2.86,
                                             ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(0, spcAxisSize, 2.86,
                                             ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);
  // lambdamax > lambda[kstart] & no igm curv -> err
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 0.5,
                                             ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);

  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tplStar.m_computedMeiksingCoeff.size() == spcAxisSize);
  BOOST_CHECK(tplStar.m_computedDustCoeff.size() == spcAxisSize);

  // InitIsmIgmConfig with lambdarange & redshift -> range outside spectral axis
  TFloat64Range lbdaRange(1, 860);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(lbdaRange, 2.86,
                                             ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);

  lbdaRange.SetBegin(spcAxisList[0]);
  lbdaRange.SetEnd(spcAxisList[spcAxisSize - 1]);
  tplStar.InitIsmIgmConfig(lbdaRange, 2.86, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tplStar.m_computedMeiksingCoeff.size() == spcAxisSize);
  BOOST_CHECK(tplStar.m_computedDustCoeff.size() == spcAxisSize);

  tplStar.ScaleFluxAxis(1.);
  tplStar.InitIsmIgmConfig(lbdaRange, 2.86, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  // InitIsmIgmConfig with redshift
  tplStar.InitIsmIgmConfig(2.86, ismCorrectionCalzetti, igmCorrectionMeiksin);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tplStar.m_computedMeiksingCoeff.size() == spcAxisSize);
  BOOST_CHECK(tplStar.m_computedDustCoeff.size() == spcAxisSize);
}

BOOST_AUTO_TEST_CASE(ApplyDustCoeff_test) {
  BOOST_CHECK_THROW(tplStar.ApplyDustCoeff(-1), GlobalException);
  BOOST_CHECK_THROW(tplStar.GetIsmCoeff(), GlobalException);

  tplStar.InitIsmIgmConfig(1, 2, 2.86, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);

  bool res = tplStar.ApplyDustCoeff(-1);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  res = tplStar.ApplyDustCoeff(-2);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  res = tplStar.ApplyDustCoeff(2);
  Int32 coeff1 =
      (Int32)(2 * ismCorrectionCalzetti->m_dataCalzetti.size() +
              round(spcAxisList[1] - ismCorrectionCalzetti->getLambdaMin()));
  Int32 coeff2 =
      (Int32)(2 * ismCorrectionCalzetti->m_dataCalzetti.size() +
              round(spcAxisList[2] - ismCorrectionCalzetti->getLambdaMin()));
  BOOST_CHECK(tplStar.GetIsmCoeff() == 2);
  BOOST_CHECK_CLOSE(tplStar.GetcomputedDustCoeffs()[1],
                    ismCorrectionCalzetti->m_dataDustCoeff[coeff1], 1e-12);
  BOOST_CHECK_CLOSE(tplStar.GetcomputedDustCoeffs()[2],
                    ismCorrectionCalzetti->m_dataDustCoeff[coeff2], 1e-12);
}

BOOST_AUTO_TEST_CASE(ApplyMeiksinCoeff_test) {
  BOOST_CHECK_THROW(tplStar.ApplyMeiksinCoeff(-1), GlobalException);

  BOOST_CHECK_THROW(tplStar.GetIgmCoeff(), GlobalException);
  BOOST_CHECK_THROW(tplStar.GetIgmEndIndex(), GlobalException);
  Int32 begin;
  Int32 ismEnd;
  BOOST_CHECK_THROW(tplStar.GetIsmIgmRangeIndex(begin, ismEnd),
                    GlobalException);

  tplStar.InitIsmIgmConfig(spcAxisSize - 2, spcAxisSize - 1, 2.86,
                           ismCorrectionCalzetti, igmCorrectionMeiksin);

  // m_Igm_kend == -1
  bool res = tplStar.ApplyMeiksinCoeff(-1);
  BOOST_CHECK(res == false);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);
  res = tplStar.ApplyMeiksinCoeff(1);
  BOOST_CHECK(res == false);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);

  // m_Igm_kend == 0
  res = tplStar.ApplyMeiksinCoeff(-1);
  BOOST_CHECK(res == true);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  // m_meiksinIdx == -2
  TAxisSampleList::const_iterator it =
      std::upper_bound(spcAxisList.begin(), spcAxisList.end(), 1216.03);
  Int32 lastIdx = it - 1 - spcAxisList.begin();
  res = tplStar.ApplyMeiksinCoeff(-2);

  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tplStar.GetIgmCoeff() == -2);
  BOOST_CHECK(tplStar.GetIgmEndIndex() == lastIdx);
  tplStar.GetIsmIgmRangeIndex(begin, ismEnd);
  BOOST_CHECK(begin == 0);
  BOOST_CHECK(ismEnd == spcAxisSize - 1);

  // test OK
  res = tplStar.ApplyMeiksinCoeff(1);
  BOOST_CHECK(res == true);
  TFloat64List meiskinCoeff = tplStar.GetcomputedMeiksinCoeffs();
  for (size_t i = 0; i < meiskinCoeff.size(); i++)
    BOOST_CHECK(meiskinCoeff[i] * fluxAxisList[i] == tplStar.GetFluxAxis()[i]);
}

BOOST_AUTO_TEST_CASE(Getter_Setter_test) {
  tplStar.SetSpectralAndFluxAxes(spcAxis, fluxAxis);
  tplStar.SetContinuumEstimationMethod("zero");

  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);

  BOOST_CHECK(tplStar.CheckIsmIgmEnabled() == true);
  tplStar.SetType(CSpectrum::EType::nType_noContinuum);
  BOOST_CHECK(tplStar.CheckIsmIgmEnabled() == false);

  CSpectrumSpectralAxis spectralAxis2(spcAxisList);
  tplStar.SetSpectralAxis(spectralAxis2);
  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  const CTemplate tpl2 = tplStar;
  BOOST_CHECK_THROW(tpl2.SetType(CSpectrum::EType::nType_continuumOnly),
                    GlobalException);
  CSpectrumFluxAxis fluxAxis2(fluxAxisList);
  tplStar.SetFluxAxis(fluxAxis2);
  tplStar.SetSpectralAxis(std::move(spcAxisList));

  // GetIgmEndIndex
  CTemplate tpl3("name", "category", tplStar.GetSpectralAxis(), fluxAxisList);
  BOOST_CHECK_THROW(tpl3.GetIgmEndIndex(0, 2), GlobalException);

  tpl3.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86, ismCorrectionCalzetti,
                        igmCorrectionMeiksin);

  TFloat64List spcAxisList = spectralAxis2.GetSamplesVector();
  TAxisSampleList::const_iterator it =
      std::upper_bound(spcAxisList.begin(), spcAxisList.end(), 1216.03);
  Int32 lastIdx = it - 1 - spcAxisList.begin();
  Int32 igmEndIndex = tpl3.GetIgmEndIndex(0, spcAxisSize - 1);
  BOOST_CHECK(igmEndIndex == lastIdx);
  igmEndIndex = tpl3.GetIgmEndIndex(0, 2);
  BOOST_CHECK(igmEndIndex == 2);
  igmEndIndex = tpl3.GetIgmEndIndex(0, 1);
  BOOST_CHECK(igmEndIndex == 1);

  tpl3.InitIsmIgmConfig(spcAxisSize - 2, spcAxisSize - 1, 2.86,
                        ismCorrectionCalzetti, igmCorrectionMeiksin);
  igmEndIndex = tpl3.GetIgmEndIndex(spcAxisSize - 2, spcAxisSize - 1);
  BOOST_CHECK(igmEndIndex == -1);

  // GetIsmIdxList
  CTemplate tpl4("name", "category", tplStar.GetSpectralAxis(), fluxAxisList);
  TInt32List ebmvList;
  BOOST_CHECK_THROW(tpl4.GetIsmIdxList(1, 1), GlobalException);

  tpl4.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86, ismCorrectionCalzetti,
                        igmCorrectionMeiksin);

  ebmvList = tpl4.GetIsmIdxList(0, 1);
  BOOST_CHECK(ebmvList.size() == 1);
  BOOST_CHECK(ebmvList[0] == -1);

  ebmvList = tpl4.GetIsmIdxList(1, -1);
  BOOST_CHECK(ebmvList.size() == 10);
  TInt32List ref_list = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  BOOST_CHECK(ebmvList == ref_list);

  ebmvList = tpl4.GetIsmIdxList(1, 3);
  BOOST_CHECK(ebmvList.size() == 1);
  BOOST_CHECK(ebmvList[0] == 3);

  // GetIgmIdxList
  CTemplate tpl5("name", "category", tplStar.GetSpectralAxis(), fluxAxisList);
  TInt32List meiksinList;
  BOOST_CHECK_THROW(tpl5.GetIgmIdxList(1, 1), GlobalException);

  tpl5.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86, ismCorrectionCalzetti,
                        igmCorrectionMeiksin);

  meiksinList = tpl5.GetIgmIdxList(0, 1);
  BOOST_CHECK(meiksinList.size() == 1);
  BOOST_CHECK(meiksinList[0] == -1);

  meiksinList = tpl5.GetIgmIdxList(1, -1);
  BOOST_CHECK(meiksinList.size() == idxCount);
  ref_list = {0, 1};
  BOOST_CHECK(meiksinList == ref_list);

  meiksinList = tpl5.GetIgmIdxList(1, 3);
  BOOST_CHECK(meiksinList.size() == 1);
  BOOST_CHECK(meiksinList[0] == 3);

  // GetIsmIgmIdxList
  TInt32List ebmvList2;
  TInt32List meiksinList2;
  tpl5.GetIsmIgmIdxList(1, 1, meiksinList2, ebmvList2, 3, 3);
  BOOST_CHECK(meiksinList2 == meiksinList);
  BOOST_CHECK(ebmvList2 == ebmvList);
}

BOOST_AUTO_TEST_SUITE_END()

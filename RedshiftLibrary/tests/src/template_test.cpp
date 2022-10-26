
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
#include "test-config.h"

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include <cmath>

using namespace NSEpic;

TFloat64List spectralList = {1213, 1214, 1215, 1216, 1217, 1218};
TFloat64List fluxList = {1e-2, 1.5e-2, 1e-2, 1.2e-2, 1e-2, 1.1e-2};
TFloat64List noiseList = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};
TFloat64List maskList = {0, 1, 1, 1, 1, 1};

class fixture_TemplateTest {
public:
  CTemplate tplStar = fixture_TemplateStarLight().tplStarLight;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      fixture_MeiskinCorrection().igmCorrectionMeiksin;
  CSpectrumSpectralAxis spcAxis = fixture_SpectralAxisLight().spcAxisLight;
  CSpectrumFluxAxis fluxAxis = fixture_FluxAxisLight().fluxAxisLight;
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

  BOOST_CHECK(tplStar.GetSampleCount() == mySpectralListLight.size());
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  // copy and copy assignement
  CTemplate tpl4(tplStar);
  BOOST_CHECK(tpl4.GetSampleCount() == mySpectralListLight.size());
  BOOST_CHECK(tpl4.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  CTemplate tpl5(tplStar);
  tpl5.GetRawFluxAxis_().GetSamplesVector().pop_back();
  BOOST_CHECK_THROW(CTemplate tpl6(tpl5), GlobalException);

  CTemplate tpl7;
  tpl7 = tplStar;
  BOOST_CHECK(tpl7.GetSampleCount() == mySpectralListLight.size());
  BOOST_CHECK(tpl7.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  // move and move assignement
  CTemplate tpl8(std::move(tplStar));
  BOOST_CHECK(tpl8.GetSampleCount() == mySpectralListLight.size());
  BOOST_CHECK(tpl8.GetFluxAxis().GetSamplesVector() == myFluxListLight);
  BOOST_CHECK(tplStar.GetSampleCount() == 0);

  BOOST_CHECK_THROW(CTemplate tpl8b(std::move(tpl5)), GlobalException);

  CTemplate tpl9;
  tpl9 = std::move(tpl8);
  BOOST_CHECK(tpl9.GetSampleCount() == mySpectralListLight.size());
  BOOST_CHECK(tpl9.GetFluxAxis().GetSamplesVector() == myFluxListLight);
  BOOST_CHECK(tpl8.GetSampleCount() == 0);

  // InitIsmIgmConfig
  tpl9.InitIsmIgmConfig(0, 2, 1214, ismCorrectionCalzetti,
                        igmCorrectionMeiksin);
  CTemplate tpl10(tpl9, maskList);
  BOOST_CHECK(tpl10.m_NoIsmIgmFluxAxis.GetSamplesCount() == 5);
  BOOST_CHECK(tpl10.m_NoIsmIgmFluxAxis[0] == 1.5e-2);
  BOOST_CHECK(tpl10.m_computedDustCoeff.size() == 5);
  BOOST_CHECK(tpl10.m_computedDustCoeff[0] == 1.);
  BOOST_CHECK(tpl10.m_computedMeiksingCoeff.size() == 5);
  BOOST_CHECK(tpl10.m_computedMeiksingCoeff[0] == 1.);

  TFloat64List maskList2 = {0, 0, 0, 1, 1, 0};
  CTemplate tpl11(tpl9, maskList2);
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
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(1, 2, 1214), GlobalException);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(8, 12, 1214, ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(0, 8, 1214, ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(0, 2, 0.5, ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);

  tplStar.InitIsmIgmConfig(1, 2, 1214, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == myFluxListLight);
  BOOST_CHECK(tplStar.m_computedMeiksingCoeff.size() ==
              mySpectralListLight.size());
  BOOST_CHECK(tplStar.m_computedDustCoeff.size() == mySpectralListLight.size());

  // InitIsmIgmConfig with lambdarange & redshift
  //   CTemplate tpl2 = ctx.createTemplate();
  TFloat64Range lbdaRange(7, 10);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(lbdaRange, 1214,
                                             ismCorrectionCalzetti,
                                             igmCorrectionMeiksin),
                    GlobalException);

  lbdaRange.SetBegin(1214);
  lbdaRange.SetEnd(1218);
  tplStar.InitIsmIgmConfig(lbdaRange, 1214, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == myFluxListLight);
  BOOST_CHECK(tplStar.m_computedMeiksingCoeff.size() ==
              mySpectralListLight.size());
  BOOST_CHECK(tplStar.m_computedDustCoeff.size() == mySpectralListLight.size());

  tplStar.ScaleFluxAxis(2);
  tplStar.InitIsmIgmConfig(lbdaRange, 2428, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  TFloat64List fluxList2 = {2e-2, 4e-2, 6e-2, 8e-2, 10e-2, 12e-2};
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxList2);

  // InitIsmIgmConfig with redshift
  //   CTemplate tpl3 = ctx.createTemplate();
  tplStar.ScaleFluxAxis(0.5);
  tplStar.InitIsmIgmConfig(2428, ismCorrectionCalzetti, igmCorrectionMeiksin);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == myFluxListLight);
  BOOST_CHECK(tplStar.m_computedMeiksingCoeff.size() ==
              mySpectralListLight.size());
  BOOST_CHECK(tplStar.m_computedDustCoeff.size() == mySpectralListLight.size());
}

BOOST_AUTO_TEST_CASE(ApplyDustCoeff_test) {
  BOOST_CHECK_THROW(tplStar.ApplyDustCoeff(-1), GlobalException);
  BOOST_CHECK_THROW(tplStar.GetIsmCoeff(), GlobalException);

  tplStar.InitIsmIgmConfig(1, 2, 1214, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);

  bool res = tplStar.ApplyDustCoeff(-1);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  res = tplStar.ApplyDustCoeff(-2);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  res = tplStar.ApplyDustCoeff(2);
  Int32 coeff1 =
      (Int32)(2 * ismCorrectionCalzetti->m_dataCalzetti.size() +
              mySpectralListLight[1] - ismCorrectionCalzetti->getLambdaMin());
  Int32 coeff2 =
      (Int32)(2 * ismCorrectionCalzetti->m_dataCalzetti.size() +
              mySpectralListLight[2] - ismCorrectionCalzetti->getLambdaMin());
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

  CSpectrumSpectralAxis spectralAxis2({1213, 1214, 1215, 1216, 1217, 1218});
  tplStar.SetSpectralAxis(spectralAxis2);
  tplStar.InitIsmIgmConfig(4, 5, 1213., ismCorrectionCalzetti,
                           igmCorrectionMeiksin);

  // m_Igm_kend == -1
  bool res = tplStar.ApplyMeiksinCoeff(-1);
  BOOST_CHECK(res == false);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == myFluxListLight);
  res = tplStar.ApplyMeiksinCoeff(1);
  BOOST_CHECK(res == false);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  tplStar.InitIsmIgmConfig(0, 2, 1214., ismCorrectionCalzetti,
                           igmCorrectionMeiksin);

  // m_Igm_kend == 2
  res = tplStar.ApplyMeiksinCoeff(-1);
  BOOST_CHECK(res == true);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  // m_meiksinIdx == -2
  res = tplStar.ApplyMeiksinCoeff(-2);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == myFluxListLight);

  BOOST_CHECK(tplStar.GetIgmCoeff() == -2);
  BOOST_CHECK(tplStar.GetIgmEndIndex() == 2);
  tplStar.GetIsmIgmRangeIndex(begin, ismEnd);
  BOOST_CHECK(begin == 0);
  BOOST_CHECK(ismEnd == 2);

  // test OK
  res = tplStar.ApplyMeiksinCoeff(1);
  BOOST_CHECK(res == true);
  TFloat64List meiskinCoeff = tplStar.GetcomputedMeiksinCoeffs();
  for (size_t i = 0; i < meiskinCoeff.size(); i++)
    BOOST_CHECK(meiskinCoeff[i] * myFluxListLight[i] ==
                tplStar.GetFluxAxis()[i]);
}

BOOST_AUTO_TEST_CASE(Getter_Setter_test) {
  tplStar.SetSpectralAndFluxAxes(spcAxis, fluxAxis);
  tplStar.SetContinuumEstimationMethod("zero");

  tplStar.InitIsmIgmConfig(0, 2, 1214, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);

  BOOST_CHECK(tplStar.CheckIsmIgmEnabled() == true);
  tplStar.SetType(CSpectrum::EType::nType_noContinuum);
  BOOST_CHECK(tplStar.CheckIsmIgmEnabled() == false);

  CSpectrumSpectralAxis spectralAxis2(mySpectralListLight);
  tplStar.SetSpectralAxis(spectralAxis2);
  tplStar.InitIsmIgmConfig(0, 2, 1214, ismCorrectionCalzetti,
                           igmCorrectionMeiksin);
  const CTemplate tpl2 = tplStar;
  BOOST_CHECK_THROW(tpl2.SetType(CSpectrum::EType::nType_continuumOnly),
                    GlobalException);
  CSpectrumFluxAxis fluxAxis2(myFluxListLight);
  tplStar.SetFluxAxis(fluxAxis2);
  tplStar.SetSpectralAxis(std::move(mySpectralListLight));

  // GetIgmEndIndex
  CTemplate tpl3("name", "category", tplStar.GetSpectralAxis(),
                 myFluxListLight);
  BOOST_CHECK_THROW(tpl3.GetIgmEndIndex(0, 2), GlobalException);

  tpl3.InitIsmIgmConfig(0, 2, 1214, ismCorrectionCalzetti,
                        igmCorrectionMeiksin);

  Int32 igmEndIndex = tpl3.GetIgmEndIndex(0, 4);
  BOOST_CHECK(igmEndIndex == 3);
  igmEndIndex = tpl3.GetIgmEndIndex(0, 2);
  BOOST_CHECK(igmEndIndex == 2);
  igmEndIndex = tpl3.GetIgmEndIndex(0, 1);
  BOOST_CHECK(igmEndIndex == 1);
  igmEndIndex = tpl3.GetIgmEndIndex(4, 5);
  BOOST_CHECK(igmEndIndex == -1);

  // GetIsmIdxList
  CTemplate tpl4("name", "category", tplStar.GetSpectralAxis(),
                 myFluxListLight);
  TInt32List ebmvList;
  BOOST_CHECK_THROW(tpl4.GetIsmIdxList(1, 1), GlobalException);

  tpl4.InitIsmIgmConfig(0, 2, 1214, ismCorrectionCalzetti,
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
  CTemplate tpl5("name", "category", tplStar.GetSpectralAxis(),
                 myFluxListLight);
  TInt32List meiksinList;
  BOOST_CHECK_THROW(tpl5.GetIgmIdxList(1, 1), GlobalException);

  tpl5.InitIsmIgmConfig(0, 2, 1214, ismCorrectionCalzetti,
                        igmCorrectionMeiksin);

  meiksinList = tpl5.GetIgmIdxList(0, 1);
  BOOST_CHECK(meiksinList.size() == 1);
  BOOST_CHECK(meiksinList[0] == -1);

  meiksinList = tpl5.GetIgmIdxList(1, -1);
  BOOST_CHECK(meiksinList.size() == 7);
  ref_list = {0, 1, 2, 3, 4, 5, 6};
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


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
#include <cmath>

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "tests/src/tool/inputContextLight.h"

using namespace NSEpic;

class fixture_TemplateTest {
public:
  fixture_TemplateTest() {
    igmCorrectionMeiksin->convolveByLSF(LSF,
                                        fixture_MeiskinCorrection().lbdaRange);
  }
  std::shared_ptr<CScopeStack> scopeStack = std::make_shared<CScopeStack>();
  CTemplate tplStar = fixture_TemplateStar().tplStar;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      fixture_CalzettiCorrection().ismCorrectionCalzetti;
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

  Context.setFluxCorrectionCalzetti(ismCorrectionCalzetti);
  Context.setFluxCorrectionMeiksin(igmCorrectionMeiksin);

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
  BOOST_CHECK_THROW(CTemplate tpl6(tpl5), AmzException);

  CTemplate tpl7;
  tpl7 = tplStar;
  BOOST_CHECK(tpl7.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tpl7.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  // move and move assignement
  CTemplate tpl8(std::move(tplStar));
  BOOST_CHECK(tpl8.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tpl8.GetFluxAxis().GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tplStar.GetSampleCount() == 0);

  BOOST_CHECK_THROW(CTemplate tpl8b(std::move(tpl5)), AmzException);

  CTemplate tpl9;
  tpl9 = std::move(tpl8);
  BOOST_CHECK(tpl9.GetSampleCount() == spcAxisSize);
  BOOST_CHECK(tpl9.GetFluxAxis().GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(tpl8.GetSampleCount() == 0);

  // InitIsmIgmConfig
  TFloat64List maskList(spcAxisSize, 1);
  maskList[0] = 0;
  tpl9.InitIsmIgmConfig(0, 40, 2.8);
  CTemplate tpl10(tpl9, maskList);
  BOOST_CHECK(tpl10.m_NoIsmIgmFluxAxis.GetSamplesCount() == spcAxisSize - 1);
  BOOST_CHECK(tpl10.m_NoIsmIgmFluxAxis[0] == fluxAxisList[1]);
  BOOST_CHECK(ssize(tpl10.m_computedDustCoeff) == spcAxisSize - 1);
  BOOST_CHECK(tpl10.m_computedDustCoeff[0] == 1.);
  BOOST_CHECK(ssize(tpl10.m_computedDustCoeff) == spcAxisSize - 1);
  BOOST_CHECK(ssize(tpl10.m_computedMeiksingCoeff) == spcAxisSize - 1);
  BOOST_CHECK(tpl10.m_computedMeiksingCoeff[0] == 1.);

  TFloat64List maskList2(spcAxisSize, 0);
  tpl7.InitIsmIgmConfig(0, 2, 2.8);
  maskList2[3] = 1.;
  CTemplate tpl11(tpl7, maskList2);
  BOOST_CHECK(tpl11.m_NoIsmIgmFluxAxis.GetSamplesCount() == 0);
  BOOST_CHECK(tpl11.m_computedDustCoeff.size() == 0);
  BOOST_CHECK(tpl11.m_computedMeiksingCoeff.size() == 0);
  Context.reset();
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
  Context.setFluxCorrectionCalzetti(ismCorrectionCalzetti);
  Context.setFluxCorrectionMeiksin(igmCorrectionMeiksin);

  // InitIsmIgmConfig with kstart, kend & redshift
  BOOST_CHECK_THROW(
      tplStar.InitIsmIgmConfig(spcAxisSize, spcAxisSize + 2, 2.86),
      AmzException);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(0, spcAxisSize, 2.86),
                    AmzException);
  // lambdamax > lambda[kstart] & no igm curv -> err
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 0.5),
                    AmzException);

  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(ssize(tplStar.m_computedMeiksingCoeff) == spcAxisSize);
  BOOST_CHECK(ssize(tplStar.m_computedDustCoeff) == spcAxisSize);

  // InitIsmIgmConfig with lambdarange & redshift -> range outside spectral axis
  TFloat64Range lbdaRange(1, 860);
  BOOST_CHECK_THROW(tplStar.InitIsmIgmConfig(lbdaRange, 2.86), AmzException);

  lbdaRange.SetBegin(spcAxisList[0]);
  lbdaRange.SetEnd(spcAxisList[spcAxisSize - 1]);
  tplStar.InitIsmIgmConfig(lbdaRange, 2.86);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(ssize(tplStar.m_computedMeiksingCoeff) == spcAxisSize);
  BOOST_CHECK(ssize(tplStar.m_computedDustCoeff) == spcAxisSize);

  tplStar.ApplyAmplitude(1.);
  tplStar.InitIsmIgmConfig(lbdaRange, 2.86);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  // InitIsmIgmConfig with redshift
  tplStar.InitIsmIgmConfig(2.86);
  BOOST_CHECK(tplStar.m_NoIsmIgmFluxAxis.GetSamplesVector() == fluxAxisList);
  BOOST_CHECK(ssize(tplStar.m_computedMeiksingCoeff) == spcAxisSize);
  BOOST_CHECK(ssize(tplStar.m_computedDustCoeff) == spcAxisSize);
  Context.reset();
}

BOOST_AUTO_TEST_CASE(ApplyDustCoeff_test) {
  Context.setFluxCorrectionCalzetti(ismCorrectionCalzetti);
  Context.setFluxCorrectionMeiksin(igmCorrectionMeiksin);

  BOOST_CHECK_THROW(tplStar.ApplyDustCoeff(-1), AmzException);
  BOOST_CHECK_THROW(tplStar.GetIsmCoeff(), AmzException);

  tplStar.InitIsmIgmConfig(1, 2, 2.86);

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
  BOOST_CHECK_THROW(tplStar.ApplyMeiksinCoeff(-1), AmzException);

  BOOST_CHECK_THROW(tplStar.GetIgmCoeff(), AmzException);
  BOOST_CHECK_THROW(tplStar.GetIgmEndIndex(), AmzException);
  BOOST_CHECK_THROW(tplStar.GetIsmIgmRangeIndex(), AmzException);

  tplStar.InitIsmIgmConfig(spcAxisSize - 2, spcAxisSize - 1, 2.86);

  // m_Igm_kend == -1
  bool res = tplStar.ApplyMeiksinCoeff(-1);
  BOOST_CHECK(res == false);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);
  res = tplStar.ApplyMeiksinCoeff(1);
  BOOST_CHECK(res == false);
  BOOST_CHECK(tplStar.GetFluxAxis().GetSamplesVector() == fluxAxisList);

  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86);

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
  auto const &[begin, ismEnd] = tplStar.GetIsmIgmRangeIndex();
  BOOST_CHECK(begin == 0);
  BOOST_CHECK(ismEnd == spcAxisSize - 1);

  // test OK
  res = tplStar.ApplyMeiksinCoeff(1);
  BOOST_CHECK(res == true);
  TFloat64List meiskinCoeff = tplStar.GetcomputedMeiksinCoeffs();
  for (size_t i = 0; i < meiskinCoeff.size(); i++)
    BOOST_CHECK(meiskinCoeff[i] * fluxAxisList[i] == tplStar.GetFluxAxis()[i]);

  Context.reset();
}

BOOST_AUTO_TEST_CASE(Getter_Setter_test) {
  Context.setFluxCorrectionCalzetti(ismCorrectionCalzetti);
  Context.setFluxCorrectionMeiksin(igmCorrectionMeiksin);

  tplStar.SetSpectralAndFluxAxes(spcAxis, fluxAxis);
  tplStar.SetContinuumEstimationMethod("zero");

  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86);

  BOOST_CHECK(tplStar.CheckIsmIgmEnabled() == true);
  tplStar.SetType(CSpectrum::EType::noContinuum);
  BOOST_CHECK(tplStar.CheckIsmIgmEnabled() == false);

  CSpectrumSpectralAxis spectralAxis2(spcAxisList);
  tplStar.SetSpectralAxis(spectralAxis2);
  tplStar.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86);
  const CTemplate tpl2 = tplStar;
  BOOST_CHECK_THROW(tpl2.SetType(CSpectrum::EType::continuumOnly),
                    AmzException);
  CSpectrumFluxAxis fluxAxis2(fluxAxisList);
  tplStar.SetFluxAxis(fluxAxis2);
  tplStar.SetSpectralAxis(std::move(spcAxisList));

  // GetIgmEndIndex
  CTemplate tpl3("name", "category", tplStar.GetSpectralAxis(), fluxAxisList);
  BOOST_CHECK_THROW(tpl3.GetIgmEndIndex(0, 2), AmzException);

  tpl3.InitIsmIgmConfig(0, spcAxisSize - 1, 2.86);

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

  tpl3.InitIsmIgmConfig(spcAxisSize - 2, spcAxisSize - 1, 2.86);
  igmEndIndex = tpl3.GetIgmEndIndex(spcAxisSize - 2, spcAxisSize - 1);
  BOOST_CHECK(igmEndIndex == -1);

  // GetIsmIgmIdxList
  TIgmIsmIdxs igmIsmIdxs = tpl3.GetIsmIgmIdxList(1, 1, 3, 3);
  BOOST_CHECK(igmIsmIdxs.igmIdxs == TInt32List{3});
  BOOST_CHECK(igmIsmIdxs.ismIdxs == TInt32List{3});

  Context.reset();
}

BOOST_AUTO_TEST_SUITE_END()

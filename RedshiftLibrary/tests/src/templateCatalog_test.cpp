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

BOOST_AUTO_TEST_SUITE(TemplateCatalog)

std::shared_ptr<CTemplate> CreateTemplate(std::string name,
                                          std::string category) {
  TFloat64List array = {0., 2., 3., 6.};
  CSpectrumSpectralAxis spectralAxis(array);
  CSpectrumFluxAxis fluxAxis(array);
  return std::make_shared<CTemplate>(name, category, spectralAxis, fluxAxis);
}

CTemplateCatalog CreateCatalog(TStringList categories,
                               TStringList tplNames_cat1,
                               TStringList tplNames_cat2) {
  if (categories.size() != 2)
    throw std::runtime_error("categories len should be two");

  CTemplateCatalog catalog(0);
  for (std::string name : tplNames_cat1)
    catalog.Add(CreateTemplate(name, categories[0]));

  for (std::string name : tplNames_cat2)
    catalog.Add(CreateTemplate(name, categories[1]));
  return catalog;
}

CTemplateCatalog CreateCatalog() {
  return CreateCatalog({"galaxy", "qso"}, {"T1", "T2"}, {"T3"});
}

void AddToCatalog(CTemplateCatalog &catalog, TStringList categories,
                  TStringList tplNames_cat1, TStringList tplNames_cat2) {
  if (categories.size() != 2)
    throw std::runtime_error("categories len should be two");

  for (std::string name : tplNames_cat1)
    catalog.Add(CreateTemplate(name, categories[0]));

  for (std::string name : tplNames_cat2)
    catalog.Add(CreateTemplate(name, categories[1]));
  return;
}

// we consider that Rebinning happens first
CTemplateCatalog CreateCatalog(bool logSampling, bool ortho) {
  CTemplateCatalog catalog = CreateCatalog();
  TStringList categories = {"galaxy", "qso"};
  if (logSampling) {
    catalog.m_logsampling = 1;
    catalog.m_orthogonal = 0;
    AddToCatalog(catalog, categories, {"rebinnedT1", "rebinnedT2"},
                 {"rebinnedT3"});
  }

  if (ortho) {
    // adding orthog templates corresponding to original templates
    catalog.m_logsampling = 0;
    catalog.m_orthogonal = 1;
    AddToCatalog(catalog, categories, {"orthoT1", "orthoT2"}, {"orthoT3"});

    if (logSampling) {
      catalog.m_logsampling = 1;
      catalog.m_orthogonal = 1;
      AddToCatalog(catalog, categories, {"rebinOrthoT1", "rebinOrthoT1"},
                   {"rebinOrthoT3"});
    }
  }
  // reinit state
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  return catalog;
}

BOOST_AUTO_TEST_CASE(Add) {

  CTemplateCatalog catalog = CreateCatalog();
  BOOST_CHECK(catalog.m_logsampling == 0);
  BOOST_CHECK(catalog.m_orthogonal == 0);

  // check size
  std::string category = "galaxy";
  BOOST_CHECK(catalog.GetTemplateCount(category) == 2);

  category = "qso";
  catalog.GetTemplateCount(category);
  BOOST_CHECK(catalog.GetTemplateCount(category) == 1);

  // check throw
  category = "";
  BOOST_CHECK_THROW(catalog.Add(CreateTemplate("T4", category)),
                    GlobalException);
}

BOOST_AUTO_TEST_CASE(SetTemplate) {
  CTemplateCatalog catalog = CreateCatalog();

  // add template to galaxy
  std::string category = "galaxy";
  catalog.SetTemplate(CreateTemplate("modifiedT1", category), 0);
  BOOST_CHECK(catalog.GetTemplateCount(category) == 2);
  // check presence of newly changed template

  std::shared_ptr<const CTemplate> retreivedTemplate =
      catalog.GetTemplate(category, 0);
  BOOST_CHECK(retreivedTemplate->GetName() == "modifiedT1");
  // case out-of-index
  BOOST_CHECK_THROW(
      catalog.SetTemplate(CreateTemplate("modifiedT1", category), 3),
      std::exception);
  // negative Index
  BOOST_CHECK_THROW(
      catalog.SetTemplate(CreateTemplate("modifiedT1", category), -2),
      std::exception);

  category = "qso";
  catalog.SetTemplate(CreateTemplate("modifiedT3", category), 0);
  BOOST_CHECK(catalog.GetTemplateCount(category) == 1);
  std::shared_ptr<const CTemplate> retreivedTemplateQSO =
      catalog.GetTemplate(category, 0);
  BOOST_CHECK(retreivedTemplateQSO->GetName() == "modifiedT3");

  category = "star";
  catalog.SetTemplate(CreateTemplate("T4", category), 0);
  BOOST_CHECK(catalog.GetTemplateCount(category) == 1);
  std::shared_ptr<const CTemplate> retreivedTemplateStar =
      catalog.GetTemplate(category, 0);
  BOOST_CHECK(retreivedTemplateStar->GetName() == "T4");
}

BOOST_AUTO_TEST_CASE(Add_rebinnedTemplates) {
  CTemplateCatalog catalog = CreateCatalog();
  catalog.m_logsampling = 1;
  BOOST_CHECK(catalog.m_logsampling == 1);
  BOOST_CHECK(catalog.m_orthogonal == 0);

  // add rebinned template to galaxy
  std::string category = "galaxy";
  catalog.Add(CreateTemplate("rebinnedT1", category));
  catalog.Add(CreateTemplate("rebinnedT2", category));

  BOOST_CHECK(catalog.GetTemplateCount(category) == 2);

  category = "qso";
  BOOST_CHECK(catalog.GetTemplateCount(category) ==
              0); // no rebinned qso yet added
  catalog.Add(CreateTemplate("rebinnedT3", category));
  BOOST_CHECK(catalog.GetTemplateCount(category) ==
              1); // no rebinned qso yet added
}

BOOST_AUTO_TEST_CASE(SetTemplate_afterChangingOriginalTemplates) {
  bool logSampling = 1, ortho = 0; // no ortho in catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  // change originalTemplates --> this should empty the correspondant rebinned
  // template
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  std::string category = "galaxy";
  catalog.SetTemplate(CreateTemplate("modifiedT1", category), 0);
  category = "qso";
  catalog.SetTemplate(CreateTemplate("modifiedT3", category), 0);

  // now check if the associated rebinned templates are nullPtr
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  category = "galaxy";
  std::shared_ptr<const CTemplate> retreivedTemplate =
      catalog.GetTemplate(category, 0);
  BOOST_CHECK(retreivedTemplate == nullptr);
  BOOST_CHECK(catalog.GetTemplateCount(category) == 2); // size doesnt change

  category = "qso";
  // check purging
  BOOST_CHECK(catalog.GetNonNullTemplateCount(category) == 0);
}

BOOST_AUTO_TEST_CASE(SetTemplate_OrthogonalnRebinnedCase) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  // change originalTemplates --> this should empty the correspondant rebinned
  // template
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  std::string category = "galaxy";
  catalog.SetTemplate(CreateTemplate("modifiedT1", category), 0);
  category = "qso";
  catalog.SetTemplate(CreateTemplate("modifiedT3", category), 0);

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 1); // decreased by 1
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") ==
              0); // purged cause initial size = 1

  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 1); // purged
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);    // purged

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 1);
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
}

BOOST_AUTO_TEST_CASE(SetTemplate_OrthogonalnRebinnedCase2) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  std::string category = "galaxy";
  catalog.SetTemplate(CreateTemplate("modifiedT1", category), 0);
  category = "qso";
  catalog.SetTemplate(CreateTemplate("modifiedT3", category), 0);

  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // same size
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // same size

  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact

  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 1); // nullptr
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);    // purged
}

BOOST_AUTO_TEST_CASE(SetTemplate_OrthogonalnRebinnedCase4) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  std::string category = "galaxy";
  catalog.SetTemplate(CreateTemplate("modifiedT1", category), 0);
  category = "qso";
  catalog.SetTemplate(CreateTemplate("modifiedT3", category), 0);

  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // same size
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // same size

  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact
}

BOOST_AUTO_TEST_CASE(SetTemplate_OrthogonalnRebinnedCase3) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  std::string category = "galaxy";
  catalog.SetTemplate(CreateTemplate("modifiedT1", category), 0);
  category = "qso";
  catalog.SetTemplate(CreateTemplate("modifiedT3", category), 0);

  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // same size
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // same size

  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact

  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("galaxy") == 2); // intact
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);    // intact
}

BOOST_AUTO_TEST_CASE(Getter_test) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  // get all categories
  TStringList categories = catalog.GetCategoryList();
  BOOST_CHECK(categories[0] == "galaxy");
  BOOST_CHECK(categories[1] == "qso");
  // create a category which  doesn't exist in catalog
  categories.push_back("star");

  // Get template by name
  std::shared_ptr<const CTemplate> tpl =
      catalog.GetTemplateByName(categories, "T1");
  BOOST_CHECK(tpl->GetCategory() == "galaxy");
  BOOST_CHECK_THROW(catalog.GetTemplateByName(categories, "T4"),
                    GlobalException);

  // get all template for actual m_logsampling & m_orthogonal
  TTemplateRefList refList = catalog.GetTemplateList(categories);
  BOOST_CHECK(refList.size() == 3);
  BOOST_CHECK(refList[0]->GetName() == "T1");
  BOOST_CHECK(refList[1]->GetName() == "T2");
  BOOST_CHECK(refList[2]->GetName() == "T3");

  const CTemplateCatalog constCatalog = CreateCatalog(logSampling, ortho);
  TTemplateConstRefList constRefList2 =
      constCatalog.GetTemplateList(categories);
  BOOST_CHECK(refList.size() == 3);
  BOOST_CHECK(refList[0]->GetName() == "T1");
  BOOST_CHECK(refList[1]->GetName() == "T2");
  BOOST_CHECK(refList[2]->GetName() == "T3");

  TTemplateConstRefList constRefList =
      constCatalog.const_TTemplateRefList_cast(refList);
  BOOST_CHECK(constRefList.size() == 3);
}

BOOST_AUTO_TEST_CASE(ClearTemplate_test) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  TStringList categories = catalog.GetCategoryList();
  // create a category which  doesn't exist in catalog
  categories.push_back("star");

  // no template to clear
  catalog.ClearTemplates("star", 0, 0, 0,
                         false); // no template to clear

  // clear all templates
  catalog.ClearTemplateList("qso", 0, 0);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);

  // reset m_ortho_LSFWidth (existing category)
  CTemplateCatalog catalog2 = CreateCatalog(logSampling, ortho);
  catalog2.m_ortho_LSFWidth = 1.0;
  catalog2.ClearTemplateList("qso", 1, 0);
  catalog2.ClearTemplateList("qso", 1, 1);
  BOOST_CHECK(std::isnan(catalog2.m_ortho_LSFWidth));

  // reset m_ortho_LSFWidth (empty category)
  logSampling = 0;
  CTemplateCatalog catalog3 = CreateCatalog(logSampling, ortho);
  catalog3.m_ortho_LSFWidth = 1.0;
  catalog3.ClearTemplateList("qso", 1, 0);
  BOOST_CHECK(std::isnan(catalog3.m_ortho_LSFWidth));
}

// case 1 : clear [ortho 1, log 1] -> other unchanged
BOOST_AUTO_TEST_CASE(ClearOneTemplate_test_Case1) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  catalog.ClearTemplates("qso", 1, 1, 0, false);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
}

// case 2 : clear [ortho 0, log 0] -> all clear
BOOST_AUTO_TEST_CASE(ClearOneTemplate_test_Case2) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  catalog.ClearTemplates("qso", 0, 0, 0, false);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
}

// case 3 : clear [ortho 1,log 0] -> other unchanged
BOOST_AUTO_TEST_CASE(ClearOneTemplate_test_Case3) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  catalog.ClearTemplates("qso", 1, 0, 0, false);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
}

// case 4: clear [ortho 0, log 1] -> [ortho 1, log 1] clear
BOOST_AUTO_TEST_CASE(ClearOneTemplate_test_Case4) {
  bool logSampling = 1, ortho = 1; // create a full catalog
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);

  catalog.ClearTemplates("qso", 0, 1, 0, false);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 0;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
  catalog.m_logsampling = 0;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 1);
  catalog.m_logsampling = 1;
  catalog.m_orthogonal = 1;
  BOOST_CHECK(catalog.GetNonNullTemplateCount("qso") == 0);
}

TFloat64List spectralList = {1213, 1214, 1215, 1216, 1217, 1218};
TFloat64List fluxList = {1e-2, 1.5e-2, 1e-2, 1.2e-2, 1e-2, 1.1e-2};
TFloat64List noiseList = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};
TFloat64List maskList = {0, 1, 1, 1, 1, 1};

const std::string jsonString =
    "{\"galaxy\" : {\"method\" : \"TemplateFittingSolve\", "
    "\"TemplateFittingSolve\" : {\"dustfit\" : true, \"extinction\" : true}},"
    "\"qso\" : {\"method\" : \"TemplateFittingSolve\", "
    "\"TemplateFittingSolve\" : {\"dustfit\" : true, \"extinction\" : true}},"
    "\"templateCatalog\" : {\"continuumRemoval\" : "
    "{ \"medianKernelWidth\" : 40.0, "
    "\"medianEvenReflection\" : false, "
    "\"method\" : \"IrregularSamplingMedian\"}}}";

class MyInputContext {
public:
  std::shared_ptr<CParameterStore> paramStore;
  std::shared_ptr<CLSF> LSF;
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti;
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin;

  void InitContext() {
    TScopeStack scopeStack;
    paramStore = std::make_shared<CParameterStore>(scopeStack);
    paramStore->FromString(jsonString);

    // create LSF
    std::string lsfType = "GaussianConstantWidth";
    paramStore->Set("LSF.width", 1.09);
    std::shared_ptr<TLSFArguments> args =
        std::make_shared<TLSFGaussianConstantWidthArgs>(paramStore);
    LSF = LSFFactory.Create(lsfType, args);

    // create Calzetti correction
    TFloat64List lbda = {1214., 1215., 1216., 1217., 1218.};
    TFloat64List flux = {0.1, 0.2, 0.5, 0.3, 0.8};
    CalzettiCorrection calzettiCorr(lbda, flux);
    ismCorrectionCalzetti = std::make_shared<CSpectrumFluxCorrectionCalzetti>(
        calzettiCorr, 0., 0.1, 10);

    // create Meiksin correction
    std::vector<MeiksinCorrection> meiskinCorr;
    meiskinCorr.push_back(MeiksinCorrection(lbda, {flux, flux, flux}));
    meiskinCorr.push_back(MeiksinCorrection(lbda, {flux, flux, flux}));
    TFloat64List z_bins = {1, 2, 3};
    igmCorrectionMeiksin =
        std::make_shared<CSpectrumFluxCorrectionMeiksin>(meiskinCorr, z_bins);
    TFloat64Range lbdaRange(1214, 1216);
    igmCorrectionMeiksin->convolveByLSF(LSF, lbdaRange);
  }

  std::shared_ptr<CParameterStore> GetParameterStore() { return paramStore; }
  std::shared_ptr<CLSF> GetLSF() { return LSF; }
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> GetMeiskinCorr() {
    return igmCorrectionMeiksin;
  }
  const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
  GetCalzettiCorr() {
    return ismCorrectionCalzetti;
  }
};

BOOST_AUTO_TEST_CASE(InitContinuumRemoval_test) {
  MyInputContext ctx;
  ctx.InitContext();

  bool logSampling = 0, ortho = 0;
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);
  std::shared_ptr<CParameterStore> paramStore = ctx.GetParameterStore();
  catalog.InitContinuumRemoval(paramStore);

  TStringList categories = catalog.GetCategoryList();
  TTemplateRefList tplList = catalog.GetTemplateList(categories);
  for (auto tpl : tplList) {
    BOOST_CHECK(tpl->GetContinuumEstimationMethod() ==
                "IrregularSamplingMedian");
    BOOST_CHECK(tpl->GetMedianWinsize() == 40.);
    BOOST_CHECK(tpl->GetMedianEvenReflection() == false);
  }
}

BOOST_AUTO_TEST_CASE(SetIsmIgmCorrection_test) {
  MyInputContext ctx;
  ctx.InitContext();

  bool logSampling = 0, ortho = 0;
  CTemplateCatalog catalog = CreateCatalog(logSampling, ortho);
  std::shared_ptr<CParameterStore> paramStore = ctx.GetParameterStore();
  const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
      ismCorrectionCalzetti = ctx.GetCalzettiCorr();
  const std::shared_ptr<const CSpectrumFluxCorrectionMeiksin>
      igmCorrectionMeiksin = ctx.GetMeiskinCorr();

  TStringList categories = catalog.GetCategoryList();
  TTemplateRefList tplList = catalog.GetTemplateList(categories);
  for (auto tpl : tplList) {
    BOOST_CHECK(tpl->CalzettiInitFailed() == true);
    BOOST_CHECK(tpl->MeiksinInitFailed() == true);
  }
  catalog.SetIsmIgmCorrection(paramStore, igmCorrectionMeiksin,
                              ismCorrectionCalzetti);
  for (auto tpl : tplList) {
    BOOST_CHECK(tpl->CalzettiInitFailed() == false);
    BOOST_CHECK(tpl->MeiksinInitFailed() == false);
  }
}

BOOST_AUTO_TEST_SUITE_END()

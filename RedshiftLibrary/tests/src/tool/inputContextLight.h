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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"
#include "RedshiftLibrary/photometry/photometricband.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/processflow/scopestore.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/rebin/rebin.h"
#include "RedshiftLibrary/spectrum/rebin/rebinLinear.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "tests/src/tool/11095-58439-0081.h"
#include "tests/src/tool/230486_spectra.h"
#include "tests/src/tool/BC03_sdss_tremonti21.h"
#include "tests/src/tool/ComposantesPCA4.h"
#include "tests/src/tool/Meiksin_Var_curves_2.5.h"
#include "tests/src/tool/Meiksin_Var_curves_3.0.h"
#include "tests/src/tool/SB_calzetti_dl1.h"
#include "tests/src/tool/linecatalogamazedvacuum_H0.h"
#include "tests/src/tool/stars_templates_munari-lowt_20170105.h"

using namespace NSEpic;
// #define fixture_Context (CProcessFlowContext::GetInstance())

static TFloat64List mySpectralListLight = {1213, 1214, 1215, 1216, 1217, 1218};
static TFloat64List myFluxListLight = {1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2};
static TFloat64List myNoiseListLight = {1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4};

// Creation of paramStore
//-----------------------

class fixture_ParamStore {
public:
  fixture_ParamStore(std::string jsonString, TScopeStack &scopeStack) {
    paramStore = std::make_shared<CParameterStore>(scopeStack);
    paramStore->FromString(jsonString);
  }
  std::shared_ptr<CParameterStore> paramStore;
};

// Creation of LSFs
//-----------------

class fixture_LSFGaussianVariableWidth {
public:
  TFloat64List wave = {1213, 1214, 1215, 1216, 1217};
  TFloat64List width = {1, 2, 3, 4, 5};
  std::shared_ptr<TLSFArguments> args =
      std::make_shared<TLSFGaussianVarWidthArgs>(wave, width);
  std::shared_ptr<CLSF> LSF = LSFFactory.Create("GaussianVariableWidth", args);
};

std::string jsonString_LSFConstantRes =
    "{\"LSF\" : {\"LSFType\" : \"GaussianConstantResolution\" , \"resolution\" "
    ": 4300}}";
class fixture_LSFGaussianConstantResolution {
public:
  fixture_LSFGaussianConstantResolution(TScopeStack scopeStack) {
    std::shared_ptr<TLSFArguments> args =
        std::make_shared<TLSFGaussianConstantResolutionArgs>(
            fixture_ParamStore(jsonString_LSFConstantRes, scopeStack)
                .paramStore);
    LSF = LSFFactory.Create("GaussianConstantResolution", args);
  }
  std::shared_ptr<CLSF> LSF;
};

class fixture_LSFGaussianNISPSIM2016 {
public:
  std::shared_ptr<TLSFArguments> args = std::make_shared<TLSFArguments>();
  std::shared_ptr<CLSF> LSF = LSFFactory.Create("GaussianNISPSIM2016", args);
};

std::string jsonString_LSFGaussian =
    "{\"LSF\" : {\"LSFType\" : \"GaussianNISPVSSPSF201707\" , \"sourcesize\" : "
    "0.1}}";
class fixture_LSFGaussianNISPVSSPSF201707 {
public:
  fixture_LSFGaussianNISPVSSPSF201707(TScopeStack scopeStack) {
    std::shared_ptr<TLSFArguments> args =
        std::make_shared<TLSFGaussianNISPVSSPSF201707Args>(
            fixture_ParamStore(jsonString_LSFGaussian, scopeStack).paramStore);
    LSF = LSFFactory.Create("GaussianNISPVSSPSF201707", args);
  }
  std::shared_ptr<CLSF> LSF;
};

std::string jsonString_LSFConstantWidth =
    "{\"LSF\" : {\"LSFType\" : \"GaussianConstantWidth\" , \"width\" : "
    "1.09}}";
class fixture_LSFGaussianConstantWidth {
public:
  fixture_LSFGaussianConstantWidth(
      TScopeStack scopeStack,
      std::string jsonStr = jsonString_LSFConstantWidth) {
    std::shared_ptr<TLSFArguments> args =
        std::make_shared<TLSFGaussianConstantWidthArgs>(
            fixture_ParamStore(jsonStr, scopeStack).paramStore);
    LSF = LSFFactory.Create("GaussianConstantWidth", args);
  }
  std::shared_ptr<CLSF> LSF;
};

// Creation of Spectrum
//---------------------

// create Spectral Axis

class fixture_SpectralAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumSpectralAxis spcAxis = spectrumData.mySpectralList;
  Int32 spcAxisSize = spectrumData.mySpectralList.size();
  TFloat64List spcAxisList = spectrumData.mySpectralList;
  TFloat64List linLambdaList = spectrumData.myLinLambdaList;
};

class fixture_SpectralAxisLog {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumSpectralAxis spcAxis = spectrumData.mySpectralListLog;
  Int32 spcAxisSize = spectrumData.mySpectralListLog.size();
  TFloat64List spcAxisList = spectrumData.mySpectralListLog;
};

class fixture_SpectralAxisLight {
public:
  CSpectrumSpectralAxis spcAxisLight = mySpectralListLight;
};

class fixture_SpectralAxisExtended {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumSpectralAxis spcAxis = spectrumData.myExtendedLambdaList;
  Int32 spcAxisSize = spectrumData.myExtendedLambdaList.size();
  TFloat64List spcAxisList = spectrumData.myExtendedLambdaList;
};

class fixture_SpectralAxisQso {
public:
  fixture_spectralQsoData spcQsoData;
  CSpectrumSpectralAxis spcAxis = spcQsoData.mySpectralList;
};

// create Noise Axis
class fixture_NoiseAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumNoiseAxis noiseAxis = spectrumData.myNoiseList;
  TFloat64List noiseList = spectrumData.myNoiseList;
};

class fixture_NoiseAxisLight {
public:
  CSpectrumNoiseAxis noiseAxisLight = myNoiseListLight;
};

class fixture_NoiseAxisExtended {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumNoiseAxis noiseAxis = spectrumData.myExtendedNoiseList;
};

class fixture_NoiseAxisQso {
public:
  fixture_spectralQsoData spcQsoData;
  CSpectrumSpectralAxis noiseAxis = spcQsoData.myNoiseList;
};

// create Flux Axis
class fixture_FluxAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumFluxAxis fluxAxis =
      CSpectrumFluxAxis(spectrumData.myFluxList, fixture_NoiseAxis().noiseAxis);
  TFloat64List fluxAxisList = spectrumData.myFluxList;
};

class fixture_FluxAxisLight {
public:
  CSpectrumFluxAxis fluxAxisLight = CSpectrumFluxAxis(
      myFluxListLight, fixture_NoiseAxisLight().noiseAxisLight);
};

class fixture_FluxAxisExtended {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumFluxAxis fluxAxis = CSpectrumFluxAxis(
      spectrumData.myExtendedFluxList, fixture_NoiseAxisExtended().noiseAxis);
  TFloat64List fluxAxisList = spectrumData.myExtendedFluxList;
};

class fixture_FluxAxisQso {
public:
  fixture_spectralQsoData spcQsoData;
  CSpectrumFluxAxis fluxAxis = CSpectrumFluxAxis(
      spcQsoData.myFluxList, fixture_NoiseAxisQso().noiseAxis);
};

// create Spectrum
class fixture_Spectrum {
public:
  CSpectrum spc =
      CSpectrum(fixture_SpectralAxis().spcAxis, fixture_FluxAxis().fluxAxis);
};

class fixture_SpectrumWithLSF {
public:
  CSpectrum spcWithVariableWidthLSF =
      CSpectrum(fixture_SpectralAxis().spcAxis, fixture_FluxAxis().fluxAxis,
                fixture_LSFGaussianVariableWidth().LSF);
};

class fixture_SpectrumLight {
public:
  CSpectrum spcLight = CSpectrum(fixture_SpectralAxisLight().spcAxisLight,
                                 fixture_FluxAxisLight().fluxAxisLight);
};

class fixture_SharedSpectrumExtended {
public:
  std::shared_ptr<CSpectrum> spc =
      std::make_shared<CSpectrum>(fixture_SpectralAxisExtended().spcAxis,
                                  fixture_FluxAxisExtended().fluxAxis);
};

class fixture_SharedSpectrum {
public:
  std::shared_ptr<CSpectrum> spc = std::make_shared<CSpectrum>(
      fixture_SpectralAxis().spcAxis, fixture_FluxAxis().fluxAxis);
};

class fixture_SharedSpectrumLog {
public:
  std::shared_ptr<CSpectrum> spc = std::make_shared<CSpectrum>(
      fixture_SpectralAxisLog().spcAxis, fixture_FluxAxis().fluxAxis);
};

class fixture_SharedSpectrumQso {
public:
  std::shared_ptr<CSpectrum> spc = std::make_shared<CSpectrum>(
      fixture_SpectralAxisQso().spcAxis, fixture_FluxAxisQso().fluxAxis);
};

// Creation of Calzetti corrections
//---------------------------------
class fixture_CalzettiCorrection {
public:
  fixture_CalzettiData calzettiData;
  CalzettiCorrection calzettiCorr =
      CalzettiCorrection(calzettiData.lbdaCorrCal, calzettiData.fluxCorrCal);
  std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti =
      std::make_shared<CSpectrumFluxCorrectionCalzetti>(calzettiCorr, 0., 0.1,
                                                        10);
};

// Creation of Meiskin corrections
//--------------------------------
class fixture_MeiskinCorrection {
public:
  fixture_MeiksinData25 meiksinData25;
  fixture_MeiksinData30 meiksinData30;
  std::vector<MeiksinCorrection> meiskinCorr = {
      MeiksinCorrection(meiksinData25.lbdaCorr,
                        {meiksinData25.fluxCorr1, meiksinData25.fluxCorr2}),
      MeiksinCorrection(meiksinData25.lbdaCorr,
                        {meiksinData30.fluxCorr1b, meiksinData30.fluxCorr2b})};
  TFloat64List z_bins = {2.0, 2.5, 3.0};
  Int32 idxCount = meiskinCorr[0].fluxcorr.size();
  TFloat64Range lbdaRange = TFloat64Range(4680., 4713.);
  std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin =
      std::make_shared<CSpectrumFluxCorrectionMeiksin>(meiskinCorr, z_bins);
};

// Creation of Template
//---------------------

class fixture_SharedStarTemplate {
public:
  fixture_SpectrumData spectrumData;
  std::shared_ptr<CTemplate> tpl = std::make_shared<CTemplate>(
      "star", "star", spectrumData.mySpectralList, spectrumData.myFluxList);
};

class fixture_SharedStarNotLogTemplate {
public:
  fixture_GalaxyTplData galaxyTplData;
  std::shared_ptr<CTemplate> tpl = std::make_shared<CTemplate>(
      "star", "star", galaxyTplData.myGalaxyLambdaList,
      galaxyTplData.myGalaxyFluxList);
};
class fixture_SharedGalaxyTemplate {
public:
  fixture_GalaxyTplData galaxyTplData;
  std::shared_ptr<CTemplate> tpl = std::make_shared<CTemplate>(
      "galaxy", "galaxy", galaxyTplData.myGalaxyLambdaList,
      galaxyTplData.myGalaxyFluxList);
  std::shared_ptr<CTemplate> tpl2 = std::make_shared<CTemplate>(
      "galaxy2", "galaxy", galaxyTplData.myGalaxyLambdaList2,
      galaxyTplData.myGalaxyFluxList2);
};

class fixture_SharedQsoTemplate {
public:
  fixture_tplQsoData tplQsoData;
  std::shared_ptr<CTemplate> tpl_c1 = std::make_shared<CTemplate>(
      "qso_c1", "qso", tplQsoData.lbda, tplQsoData.flux_c1);
  std::shared_ptr<CTemplate> tpl_c2 = std::make_shared<CTemplate>(
      "qso_c2", "qso", tplQsoData.lbda, tplQsoData.flux_c2);
  std::shared_ptr<CTemplate> tpl_c3 = std::make_shared<CTemplate>(
      "qso_c3", "qso", tplQsoData.lbda, tplQsoData.flux_c3);
  std::shared_ptr<CTemplate> tpl_c4 = std::make_shared<CTemplate>(
      "qso_c4", "qso", tplQsoData.lbda, tplQsoData.flux_c4);
  std::shared_ptr<CTemplate> tpl_mean = std::make_shared<CTemplate>(
      "qso_mean", "qso", tplQsoData.lbda, tplQsoData.flux_mean);
};

class fixture_TemplateStar {
public:
  fixture_GalaxyTplData galaxyTplData;
  CTemplate tplStar =
      CTemplate("tpl_star", "star", galaxyTplData.myGalaxyLambdaList,
                galaxyTplData.myGalaxyFluxList);
  Int32 spcAxisSize = galaxyTplData.myGalaxyLambdaList.size();
  TFloat64List spcAxisList = galaxyTplData.myGalaxyLambdaList;
  TFloat64List fluxAxisList = galaxyTplData.myGalaxyFluxList;
};

class fixture_TemplateGalaxy {
public:
  CTemplate tplGalaxy =
      CTemplate("tpl_galaxy", "galaxy", fixture_SpectralAxis().spcAxis,
                fixture_FluxAxis().fluxAxis);
};

// Creation of Catalog
//--------------------
class fixture_TemplateCatalog {
public:
  fixture_TemplateCatalog() {
    catalog.Add(fixture_SharedStarTemplate().tpl);
    catalog.Add(fixture_SharedGalaxyTemplate().tpl);
  }
  CTemplateCatalog catalog = CTemplateCatalog(0);
};

class fixture_sharedTemplateCatalog {
public:
  // fixture_sharedTemplateCatalog() {
  //   catalog->Add(fixture_SharedStarTemplate().tpl);
  //   catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  //   catalog->m_logsampling = 1;
  //   catalog->Add(fixture_SharedGalaxyTemplate().tpl);
  // }
  std::shared_ptr<CTemplateCatalog> catalog =
      std::make_shared<CTemplateCatalog>(0);
};

// Creation of photoBandCatalog
//-----------------------------

class fixture_PhotoBand {
public:
  fixture_PhotoBand() {
    TFloat64List trans = {.2, .5, .99, .99, .4, .1};
    TFloat64List lambda = {4680., 4690., 4695., 4700., 4710., 4712.};
    photoBand = CPhotometricBand(trans, lambda);
  }
  CPhotometricBand photoBand;
};

class fixture_PhotoBandCatalog {
public:
  fixture_PhotoBandCatalog() {
    TStringList names = {"band1", "band2"};
    photoBandCatalog->Add(names[0], fixture_PhotoBand().photoBand);
    photoBandCatalog->Add(names[1], fixture_PhotoBand().photoBand);
  }
  std::shared_ptr<CPhotBandCatalog> photoBandCatalog =
      std::make_shared<CPhotBandCatalog>();
};

class fixture_PhotoData {
public:
  fixture_PhotoData() {
    TStringList names = {"band1", "band2"};
    TFloat64List flux = {1e-14, 2e-15};
    TFloat64List fluxerr = {1e-18, 3e-18};
    photoData = std::make_shared<CPhotometricData>(names, flux, fluxerr);
  }
  std::shared_ptr<CPhotometricData> photoData;
};

// Creation of line catalog
//-------------------------

class fixture_LineCatalog {
public:
  fixture_LineCatalog() {

    for (std::size_t i = 0; i < lineCatalogData.waveLength.size(); i++) {
      lineCatalog->AddLineFromParams(
          lineCatalogData.name[i], lineCatalogData.waveLength[i],
          lineCatalogData.type[i], lineCatalogData.force[i],
          lineCatalogData.profile[i], TAsymParams(0, 0, 0),
          lineCatalogData.amplitudeGroupName[i],
          lineCatalogData.amplitudeGroupValue[i],
          lineCatalogData.dispersionVelocityGroupName[i],
          lineCatalogData.waveLengthOffset[i],
          lineCatalogData.enableFitWaveLengthOffset[i], 0,
          lineCatalogData.name[i],
          fixture_MeiskinCorrection().igmCorrectionMeiksin);
    }
  }
  fixture_LineCatalogData lineCatalogData;
  Int32 nsigmasupport = 8;
  std::shared_ptr<CLineCatalog> lineCatalog =
      std::make_shared<CLineCatalog>(nsigmasupport);
  Int32 lineCatalogSize = lineCatalogData.waveLength.size();
};

// Creation of line ratio catalog
//-------------------------------

class fixture_LineRatioTplCatalog {
public:
  std::shared_ptr<CLineCatalogsTplRatio> lineRatioTplCatalog =
      std::make_shared<CLineCatalogsTplRatio>();
};

class fixture_LineRatioCatalog {
public:
  fixture_LineRatioCatalog() {
    // 3727.09	[OII]3726	E	7.65e-15
    // 3729.88 [OII] 3729 E 7.65e-15
    lineRatioCatalog->setLineAmplitude("[OII]3729", 3729.88);
    lineRatioCatalog->setLineAmplitude("[OII]3726", 3727.09);
    lineRatioCatalog->addVelocity("em_vel", 320.0);
    lineRatioCatalog->addVelocity("abs_vel", 640.0);
    lineRatioCatalog->setAsymProfileAndParams("ASYMFIXED",
                                              TAsymParams(3.0, 1.5, 1.0));
    // TODO mettre condition enableIGM
    lineRatioCatalog->convertLineProfiles2SYMIGM(
        fixture_MeiskinCorrection().igmCorrectionMeiksin);

    lineRatioCatalog->setIsmIndex(0);
    lineRatioCatalog->setPrior(1);
  }
  std::shared_ptr<CLineRatioCatalog> lineRatioCatalog =
      std::make_shared<CLineRatioCatalog>("ratioCatalog",
                                          *fixture_LineCatalog().lineCatalog);
};

class fixture_InputContext {
public:
  fixture_InputContext(std::string jsonString,
                       std::shared_ptr<CParameterStore> paramStore) {

    std::shared_ptr<CSpectrum> spc;
    ctx = std::make_shared<CInputContext>(paramStore);

    spc = fixture_SharedSpectrumLog().spc;

    ctx->m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(
        paramStore->Get<TFloat64Range>("lambdarange")));
    ctx->addSpectrum(spc);
    ctx->m_logGridStep = spc->GetSpectralAxis().GetlogGridStep();
    ctx->addRebinSpectrum(spc);
  }
  std::shared_ptr<CInputContext> ctx;
};
class fixture_InputContext2 {
public:
  fixture_InputContext2(std::string jsonString,
                        std::shared_ptr<CParameterStore> paramStore) {
    std::shared_ptr<CSpectrum> spc;
    ctx = std::make_shared<CInputContext>(paramStore);
    ctx->m_lambdaRanges.push_back(std::make_shared<TFloat64Range>(
        paramStore->Get<TFloat64Range>("lambdarange")));

    spc = fixture_SharedSpectrum().spc;
    TFloat64List spcAxis =
        fixture_SharedSpectrum().spc->GetSpectralAxis().GetSamplesVector();
    spcAxis[1] = 4.681000E+03;
    spc->SetSpectralAxis(spcAxis);
    ctx->addSpectrum(spc);
    std::map<std::string, bool> fft_processing;
    paramStore->hasToLogRebin({"galaxy"}, fft_processing);
    ctx->m_logGridStep =
        paramStore->getMinZStepForFFTProcessing(fft_processing);
  }
  std::shared_ptr<CInputContext> ctx;
};
class fixture_Context {
public:
  void loadParameterStore(std::string jsonString) {
    Context.LoadParameterStore(jsonString);
  }

  void setCorrections(
      std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin,
      std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti) {
    Context.setfluxCorrectionCalzetti(ismCorrectionCalzetti);
    Context.setfluxCorrectionMeiksin(igmCorrectionMeiksin);
  }

  void setCatalog(std::shared_ptr<CTemplateCatalog> catalog) {
    Context.setTemplateCatalog(catalog);
  }

  void setPhotoBandCatalog(std::shared_ptr<CPhotBandCatalog> photoBandCatalog) {
    Context.setPhotBandCatalog(photoBandCatalog);
  }

  void
  setLineRatioCatalogCatalog(std::string objectType,
                             std::shared_ptr<CLineCatalogsTplRatio> catalog) {
    Context.setLineRatioCatalogCatalog(objectType, catalog);
  }

  void setLineCatalog(std::string objectType, std::string method,
                      std::shared_ptr<CLineCatalog> catalog) {
    Context.setLineCatalog(objectType, method, catalog);
  }

  void addSpectrum(std::shared_ptr<CSpectrum> spc, std::shared_ptr<CLSF> LSF) {
    spc->SetLSF(LSF);
    Context.addSpectrum(spc);
  }
  void initContext() { Context.Init(); }

  void reset() { Context.reset(); }
};

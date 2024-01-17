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
#include "tests/src/tool/134845168.h"
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

class fixture_MultiSpectralAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumSpectralAxis spcAxisA = spectrumData.lambdaAList;
  CSpectrumSpectralAxis spcAxisB = spectrumData.lambdaBList;
  Int32 spcAxisASize = spectrumData.lambdaAList.size();
  Int32 spcAxisBSize = spectrumData.lambdaBList.size();
  TFloat64List spcAxisAList = spectrumData.lambdaAList;
  TFloat64List spcAxisBList = spectrumData.lambdaBList;
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

class fixture_SpectralAxisFull {
public:
  CSpectrumSpectralAxis spcAxis = lambdaE;
  Int32 spcAxisSize = lambdaE.size();
  TFloat64List spcAxisList = lambdaE;
};

// create Noise Axis
class fixture_NoiseAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumNoiseAxis noiseAxis = spectrumData.myNoiseList;
  TFloat64List noiseList = spectrumData.myNoiseList;
};

class fixture_MultiNoiseAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumNoiseAxis noiseAAxis = spectrumData.noiseAList;
  CSpectrumNoiseAxis noiseBAxis = spectrumData.noiseBList;
  TFloat64List noiseAList = spectrumData.noiseAList;
  TFloat64List noiseBList = spectrumData.noiseBList;
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

class fixture_NoiseAxisFull {
public:
  CSpectrumNoiseAxis noiseAxis = errorE;
};

// create Flux Axis
class fixture_FluxAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumFluxAxis fluxAxis =
      CSpectrumFluxAxis(spectrumData.myFluxList, fixture_NoiseAxis().noiseAxis);
  TFloat64List fluxAxisList = spectrumData.myFluxList;
};

class fixture_MultiFluxAxis {
public:
  fixture_SpectrumData spectrumData;
  CSpectrumFluxAxis fluxAAxis = CSpectrumFluxAxis(
      spectrumData.fluxAList, fixture_MultiNoiseAxis().noiseAAxis);
  CSpectrumFluxAxis fluxBAxis = CSpectrumFluxAxis(
      spectrumData.fluxBList, fixture_MultiNoiseAxis().noiseBAxis);
  TFloat64List fluxAxisAList = spectrumData.fluxAList;
  TFloat64List fluxAxisBList = spectrumData.fluxBList;
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

class fixture_FluxAxisFull {
public:
  CSpectrumFluxAxis fluxAxis =
      CSpectrumFluxAxis(fluxE, fixture_NoiseAxisFull().noiseAxis);
  TFloat64List fluxAxisList = fluxE;
};

// create Spectrum
class fixture_Spectrum {
public:
  CSpectrum spc =
      CSpectrum(fixture_SpectralAxis().spcAxis, fixture_FluxAxis().fluxAxis);
};

class fixture_MultiSpectrum {
public:
  CSpectrum spcA = CSpectrum(fixture_MultiSpectralAxis().spcAxisA,
                             fixture_MultiFluxAxis().fluxAAxis);
  CSpectrum spcB = CSpectrum(fixture_MultiSpectralAxis().spcAxisB,
                             fixture_MultiFluxAxis().fluxBAxis);
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

class fixture_SharedMultiSpectrum {
public:
  std::shared_ptr<CSpectrum> spcA = std::make_shared<CSpectrum>(
      fixture_MultiSpectralAxis().spcAxisA, fixture_MultiFluxAxis().fluxAAxis);
  std::shared_ptr<CSpectrum> spcB = std::make_shared<CSpectrum>(
      fixture_MultiSpectralAxis().spcAxisB, fixture_MultiFluxAxis().fluxBAxis);
};

class fixture_SharedSpectrumFull {
public:
  std::shared_ptr<CSpectrum> spc = std::make_shared<CSpectrum>(
      fixture_SpectralAxisFull().spcAxis, fixture_FluxAxisFull().fluxAxis);
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
  std::shared_ptr<CTemplate> tpl3 = std::make_shared<CTemplate>(
      "galaxy3", "galaxy", galaxyTplData.myGalaxyLambdaListFull,
      galaxyTplData.myGalaxyFluxListFull);
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
          lineCatalogData.enableFitWaveLengthOffset[i], i,
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

/* class fixture_FullLineCatalog {
public:
  fixture_FullLineCatalog() {

    lineCatalog->AddLineFromParams(
        "P5A", 12821.59, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 0, "P5A_12821.59_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P6A", 10941.09, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 1, "P6A_10941.09_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P7A", 10052.13, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 2, "P7A_10052.13_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P8A", 9548.59, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 3, "P8A_9548.59_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P9A", 9231.55, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 4, "P9A_9231.55_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P10A", 9017.38, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 5, "P10A_9017.38_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P11A", 8865.22, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 6, "P11A_8865.22_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CaII_t3A", 8664.5, "A", "S", "SYM", TAsymParams(0, 0, 0), "T_Ca", 13.0,
        "Abs1", 0, false, 7, "CaII_t3A_8664.5_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CaII_t2A", 8544.42, "A", "S", "SYM", TAsymParams(0, 0, 0), "T_ca",
        17.0, "Abs1", 0, false, 8, "CaII_t2A_8544.42_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CaII_t1A", 8500.35, "A", "S", "SYM", TAsymParams(0, 0, 0), "T_ca",
        16.0, "Abs1", 0, false, 9, "CaII_t1A_8500.35_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HalphaA", 6564.61, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 10, "HalphaA_6564.61_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "NaD", 5895.6, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Abs1",
        0, false, 11, "NaD_5895.6_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "MgI5175", 5176.71, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 12, "MgI5175_5176.71_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HbetaA", 4862.72, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 13, "HbetaA_4862.72_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HgammaA", 4341.58, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 14, "HgammaA_4341.58_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "GBand", 4304.57, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 15, "GBand_4304.57_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HdeltaA", 4102.81, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 16, "HdeltaA_4102.81_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CaII_H", 3969.55, "A", "S", "SYM", TAsymParams(0, 0, 0), "A_Ca", 22.0,
        "Abs1", 0, false, 17, "CaII_H_3969.55_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CaII_K", 3934.73, "A", "S", "SYM", TAsymParams(0, 0, 0), "A_Ca", 23.0,
        "Abs1", 0, false, 18, "CaII_K_3934.73_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H8A", 3890.11, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 19, "H8A_3890.11_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H9A", 3836.43, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 20, "H9A_3836.43_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H10A", 3798.93, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 21, "H10A_3798.93_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H11A", 3771.65, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 22, "H11A_3771.65_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeI", 3582.19, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 23, "FeI_3582.19_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "MgI2852", 2853.73, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 24, "MgI2852_2853.73_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "MgII2803", 2804.29, "A", "S", "SYM", TAsymParams(0, 0, 0), "A_MgII",
        0.3054, "Abs1", 0, false, 25, "MgII2803_2804.29_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "MgII2796", 2797.11, "A", "S", "SYM", TAsymParams(0, 0, 0), "A_MgII",
        0.6123, "Abs1", 0, false, 26, "MgII2796_2797.11_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII2600", 2600.87, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 27, "FeII2600_2600.87_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII2586", 2587.35, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 28, "FeII2586_2587.35_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII2382", 2383.4, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 29, "FeII2382_2383.4_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII2374", 2375.1, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 30, "FeII2374_2375.1_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII2344", 2344.84, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 31, "FeII2344_2344.84_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII2260", 2261.39, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 32, "FeII2260_2261.39_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII2249", 2250.49, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 33, "FeII2249_2250.49_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "AlIII1862", 1863.29, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 34, "AlIII1862_1863.29_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "AlIII1854", 1855.22, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 35, "AlIII1854_1855.22_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "AlII1670", 1670.77, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 36, "AlII1670_1670.77_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII1608", 1608.44, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 37, "FeII1608_1608.44_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CI", 1560.73, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Abs1",
        0, false, 38, "CI_1560.73_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CIV1550A", 1550.79, "A", "S", "SYM", TAsymParams(0, 0, 0), "A_CIV",
        0.9, "Abs1", -150, false, 39, "CIV1550A_1550.79_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CIV1548A", 1548.19, "A", "S", "SYM", TAsymParams(0, 0, 0), "A_CIV",
        1.0, "Abs1", -150, false, 40, "CIV1548A_1548.19_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SiII1526", 1526.71, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 41, "SiII1526_1526.71_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SiIV1402", 1402.77, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 42, "SiIV1402_1402.77_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SiIV1393", 1393.74, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 43, "SiIV1393_1393.74_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "NiII", 1370.5, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 44, "NiII_1370.5_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CII", 1334.53, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 45, "CII_1334.53_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OI1302", 1302.15, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 46, "OI1302_1302.15_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SiII1260", 1260.42, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 47, "SiII1260_1260.42_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "LyAA", 1216.03, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 48, "LyAA_1216.03_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SiIII1206", 1206.5, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 49, "SiIII1206_1206.5_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SiII1193", 1193.29, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 50, "SiII1193_1193.29_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SiII1190", 1190.42, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", -150, false, 51, "SiII1190_1190.42_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CIII1176a", 1175.71, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 52, "CIII1176a_1175.71_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CIII1176b", 1176.37, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 53, "CIII1176b_1176.37_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "NII1084a", 1083.99, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 54, "NII1084a_1083.99_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "NII1084b", 1084.58, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 55, "NII1084b_1084.58_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OVI1037A", 1037.62, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 56, "OVI1037A_1037.62_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OVI1031A", 1031.93, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 57, "OVI1031A_1031.93_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "LyBA", 1025.72, "A", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 58, "LyBA_1025.72_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "LyGA", 972.53, "A", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Abs1", 0, false, 59, "LyGA_972.53_A",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P5", 12821.59, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 60, "P5_12821.59_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P6", 10941.09, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 61, "P6_10941.09_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P7", 10052.13, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 62, "P7_10052.13_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P8", 9548.59, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 63, "P8_9548.59_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P9", 9231.55, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 64, "P9_9231.55_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P10", 9017.38, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 65, "P10_9017.38_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "P11", 8865.22, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 66, "P11_8865.22_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SIII9530", 9533.2, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 67, "SIII9530_9533.2_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "SIII9068", 9071.1, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 68, "SIII9068_9071.1_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "ArIII7751", 7753.2, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 69, "ArIII7751_7753.2_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OII7330", 7332.2, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 70, "OII7330_7332.2_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OII7319", 7321.0, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 71, "OII7319_7321.0_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "ArIII7136", 7138.73, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 72, "ArIII7136_7138.73_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[SII]6731", 6732.68, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 73, "[SII]6731_6732.68_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[SII]6716", 6718.29, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 74, "[SII]6716_6718.29_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[NII](doublet-1)", 6585.27, "E", "W", "SYM", TAsymParams(0, 0, 0),
        "E_NII", 2.95, "Em1", 0, false, 75, "[NII](doublet-1)_6585.27_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "Halpha", 6564.61, "E", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 76, "Halpha_6564.61_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[NII](doublet-1/2.95)", 6549.86, "E", "W", "SYM", TAsymParams(0, 0, 0),
        "E_NII", 1.0, "Em1", 0, false, 77, "[NII](doublet-1/2.95)_6549.86_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeI6494", 6496.75, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 78, "FeI6494_6496.75_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[OI]6301", 6303.05, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 79, "[OI]6301_6303.05_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HeI5876", 5877.41, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 80, "HeI5876_5877.41_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[OIII](doublet-1)", 5008.24, "E", "S", "SYM", TAsymParams(0, 0, 0),
        "E_OIII", 3.0, "Em1", 0, false, 81, "[OIII](doublet-1)_5008.24_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[OIII](doublet-1/3)", 4960.29, "E", "S", "SYM", TAsymParams(0, 0, 0),
        "E_OIII", 1.0, "Em1", 0, false, 82, "[OIII](doublet-1/3)_4960.29_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "Hbeta", 4862.72, "E", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 83, "Hbeta_4862.72_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "Hgamma", 4341.58, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 84, "Hgamma_4341.58_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "Hdelta", 4102.81, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 85, "Hdelta_4102.81_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "Hepsilon", 3971.15, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 86, "Hepsilon_3971.15_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H8", 3890.11, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 87, "H8_3890.11_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "NeIII", 3869.05, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 88, "NeIII_3869.05_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H9", 3836.43, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 89, "H9_3836.43_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H10", 3798.93, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 90, "H10_3798.93_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeII3785", 3785.47, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 91, "FeII3785_3785.47_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "H11", 3771.65, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 92, "H11_3771.65_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "FeVII", 3758.46, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 93, "FeVII_3758.46_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[OII]3729", 3729.88, "E", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 94, "[OII]3729_3729.88_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[OII]3726", 3727.09, "E", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 95, "[OII]3726_3727.09_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[NeVa]", 3426.73, "E", "W", "SYM", TAsymParams(0, 0, 0), "NeV_a", 2.76,
        "Em1", 0, false, 96, "[NeVa]_3426.73_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[NeVb]", 3346.81, "E", "W", "SYM", TAsymParams(0, 0, 0), "NeV_b", 1.0,
        "Em1", 0, false, 97, "[NeVb]_3346.81_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HeI3190", 3191.78, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 98, "HeI3190_3191.78_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "MgII", 2799.12, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 99, "MgII_2799.12_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[CIII]1907", 1906.68, "E", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 100, "[CIII]1907_1906.68_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "[CIII]1909", 1908.73, "E", "S", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 101, "[CIII]1909_1908.73_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OIII1661", 1661.26, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 102, "OIII1661_1661.26_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OIII1666", 1666.6, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 103, "OIII1666_1666.6_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HeII(doublet-3)", 1640.47, "E", "W", "SYM", TAsymParams(0, 0, 0),
        "E_HeII", 3.0, "Em1", 0, false, 104, "HeII(doublet-3)_1640.47_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "HeII(doublet-2)", 1640.33, "E", "W", "SYM", TAsymParams(0, 0, 0),
        "E_HeII", 2.0, "Em1", 0, false, 105, "HeII(doublet-2)_1640.33_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CIV1550", 1551.19, "E", "S", "SYM", TAsymParams(0, 0, 0), "E_CIV", 0.9,
        "Em1", 0, false, 106, "CIV1550_1551.19_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "CIV1548", 1548.61, "E", "S", "SYM", TAsymParams(0, 0, 0), "E_CIV", 1.0,
        "Em1", 0, false, 107, "CIV1548_1548.61_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "LyAE", 1216.03, "E", "S", "ASYMFIT", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 108, "LyAE_1216.03_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OVI1037E", 1037.62, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 109, "OVI1037E_1037.62_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "OVI1031E", 1031.93, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 110, "OVI1031E_1031.93_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "LyBE", 1025.72, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0,
        "Em1", 0, false, 111, "LyBE_1025.72_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
    lineCatalog->AddLineFromParams(
        "LyGE", 972.53, "E", "W", "SYM", TAsymParams(0, 0, 0), "-1", 1.0, "Em1",
        0, false, 112, "LyGE_972.53_E",
        fixture_MeiskinCorrection().igmCorrectionMeiksin);
  }
  Int32 nsigmasupport = 8;
  std::shared_ptr<CLineCatalog> lineCatalog =
      std::make_shared<CLineCatalog>(nsigmasupport);
  Int32 lineCatalogSize = 115;
};
 */

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
    lineRatioCatalog->setLineAmplitude(0, 3729.88);
    lineRatioCatalog->setLineAmplitude(1, 3727.09);
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

/* class fixture_LineRatioCatalogA {
public:
  fixture_LineRatioCatalogA() {

    lineRatioCatalogA->setLineAmplitude("LyGE_972.53_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("LyGA_972.53_A", 1.04067);
    lineRatioCatalogA->setLineAmplitude("LyBA_1025.72_A", 0.545213);
    lineRatioCatalogA->setLineAmplitude("LyBE_1025.72_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("OVI1031E_1031.93_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("OVI1031A_1031.93_A", 0.370479);
    lineRatioCatalogA->setLineAmplitude("OVI1037E_1037.62_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("OVI1037A_1037.62_A", 0.00660704);
    lineRatioCatalogA->setLineAmplitude("CIII1176a_1175.71_A", 0.0);
    lineRatioCatalogA->setLineAmplitude("CIII1176b_1176.37_A", 0.280219);
    lineRatioCatalogA->setLineAmplitude("SiII1190_1190.42_A", 0.092508);
    lineRatioCatalogA->setLineAmplitude("SiII1193_1193.29_A", 0.251816);
    lineRatioCatalogA->setLineAmplitude("SiIII1206_1206.5_A", 0.530322);
    lineRatioCatalogA->setLineAmplitude("LyAA_1216.03_A", 0.0);
    lineRatioCatalogA->setLineAmplitude("LyAE_1216.03_E", 3.12014e-15);
    lineRatioCatalogA->setLineAmplitude("SiII1260_1260.42_A", 0.177301);
    lineRatioCatalogA->setLineAmplitude("OI1302_1302.15_A", 0.264816);
    lineRatioCatalogA->setLineAmplitude("CII_1334.53_A", 0.201722);
    lineRatioCatalogA->setLineAmplitude("NiII_1370.5_A", 0.0118898);
    lineRatioCatalogA->setLineAmplitude("SiIV1393_1393.74_A", 0.229596);
    lineRatioCatalogA->setLineAmplitude("SiIV1402_1402.77_A", 0.136896);
    lineRatioCatalogA->setLineAmplitude("SiII1526_1526.71_A", 0.17199);
    lineRatioCatalogA->setLineAmplitude("CIV1548A_1548.19_A", 0.162787);
    lineRatioCatalogA->setLineAmplitude("CIV1548_1548.61_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("CIV1550A_1550.79_A", 0.146509);
    lineRatioCatalogA->setLineAmplitude("CIV1550_1551.19_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("CI_1560.73_A", 0.0);
    lineRatioCatalogA->setLineAmplitude("FeII1608_1608.44_A", 0.121272);
    lineRatioCatalogA->setLineAmplitude("HeII(doublet-2)_1640.33_E",
                                        1.25198e-16);
    lineRatioCatalogA->setLineAmplitude("HeII(doublet-3)_1640.47_E",
                                        1.87797e-16);
    lineRatioCatalogA->setLineAmplitude("OIII1661_1661.26_E", 1.59357e-16);
    lineRatioCatalogA->setLineAmplitude("OIII1666_1666.6_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("AlII1670_1670.77_A", 0.0796577);
    lineRatioCatalogA->setLineAmplitude("AlIII1854_1855.22_A", 0.0902249);
    lineRatioCatalogA->setLineAmplitude("AlIII1862_1863.29_A", 0.0795008);
    lineRatioCatalogA->setLineAmplitude("[CIII]1907_1906.68_E", 3.133e-16);
    lineRatioCatalogA->setLineAmplitude("[CIII]1909_1908.73_E", 3.133e-16);
    lineRatioCatalogA->setLineAmplitude("FeII2249_2250.49_A", 0.0294006);
    lineRatioCatalogA->setLineAmplitude("FeII2260_2261.39_A", 0.0628961);
    lineRatioCatalogA->setLineAmplitude("FeII2344_2344.84_A", 0.239228);
    lineRatioCatalogA->setLineAmplitude("FeII2374_2375.1_A", 0.13715);
    lineRatioCatalogA->setLineAmplitude("FeII2382_2383.4_A", 0.195162);
    lineRatioCatalogA->setLineAmplitude("FeII2586_2587.35_A", 0.189696);
    lineRatioCatalogA->setLineAmplitude("FeII2600_2600.87_A", 0.229154);
    lineRatioCatalogA->setLineAmplitude("MgII2796_2797.11_A", 0.184851);
    lineRatioCatalogA->setLineAmplitude("MgII_2799.12_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("MgII2803_2804.29_A", 0.0921989);
    lineRatioCatalogA->setLineAmplitude("MgI2852_2853.73_A", 0.0452975);
    lineRatioCatalogA->setLineAmplitude("HeI3190_3191.78_E",
                                        6.806640000000002e-18);
    lineRatioCatalogA->setLineAmplitude("[NeVb]_3346.81_E",
                                        7.053410000000001e-18);
    lineRatioCatalogA->setLineAmplitude("[NeVa]_3426.73_E", 3.10385e-17);
    lineRatioCatalogA->setLineAmplitude("FeI_3582.19_A", 0.015086);
    lineRatioCatalogA->setLineAmplitude("[OII]3726_3727.09_E", 1e-16);
    lineRatioCatalogA->setLineAmplitude("[OII]3729_3729.88_E", 1e-16);
    lineRatioCatalogA->setLineAmplitude("FeVII_3758.46_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("H11_3771.65_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("H11A_3771.65_A", 0.233965);
    lineRatioCatalogA->setLineAmplitude("FeII3785_3785.47_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("H10_3798.93_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("H10A_3798.93_A", 0.248023);
    lineRatioCatalogA->setLineAmplitude("H9A_3836.43_A", 0.281598);
    lineRatioCatalogA->setLineAmplitude("H9_3836.43_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("NeIII_3869.05_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("H8_3890.11_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("H8A_3890.11_A", 0.208883);
    lineRatioCatalogA->setLineAmplitude("CaII_K_3934.73_A", 0.155025);
    lineRatioCatalogA->setLineAmplitude("CaII_H_3969.55_A", 0.148285);
    lineRatioCatalogA->setLineAmplitude("Hepsilon_3971.15_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("HdeltaA_4102.81_A", 0.133059);
    lineRatioCatalogA->setLineAmplitude("Hdelta_4102.81_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("GBand_4304.57_A", 0.0782481);
    lineRatioCatalogA->setLineAmplitude("HgammaA_4341.58_A", 0.0287842);
    lineRatioCatalogA->setLineAmplitude("Hgamma_4341.58_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("HbetaA_4862.72_A", 0.0);
    lineRatioCatalogA->setLineAmplitude("Hbeta_4862.72_E", 4.5e-17);
    lineRatioCatalogA->setLineAmplitude("[OIII](doublet-1/3)_4960.29_E",
                                        1.94903e-16);
    lineRatioCatalogA->setLineAmplitude("[OIII](doublet-1)_5008.24_E",
                                        5.84709e-16);
    lineRatioCatalogA->setLineAmplitude("MgI5175_5176.71_A", 0.047198);
    lineRatioCatalogA->setLineAmplitude("HeI5876_5877.41_E",
                                        9.527570000000001e-18);
    lineRatioCatalogA->setLineAmplitude("NaD_5895.6_A", 0.0243732);
    lineRatioCatalogA->setLineAmplitude("[OI]6301_6303.05_E", 7.59793e-17);
    lineRatioCatalogA->setLineAmplitude("[NII](doublet-1/2.95)_6549.86_E",
                                        9.79943e-17);
    lineRatioCatalogA->setLineAmplitude("HalphaA_6564.61_A", 0.0);
    lineRatioCatalogA->setLineAmplitude("Halpha_6564.61_E", 4e-16);
    lineRatioCatalogA->setLineAmplitude("[NII](doublet-1)_6585.27_E",
                                        2.89083e-16);
    lineRatioCatalogA->setLineAmplitude("[SII]6716_6718.29_E", 3.31454e-16);
    lineRatioCatalogA->setLineAmplitude("[SII]6731_6732.68_E", 2.56641e-16);
    lineRatioCatalogA->setLineAmplitude("ArIII7136_7138.73_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("OII7319_7321.0_E", 0.0);
    lineRatioCatalogA->setLineAmplitude("OII7330_7332.2_E", 4.25514e-17);
    lineRatioCatalogA->setLineAmplitude("ArIII7751_7753.2_E", 2.72805e-17);
    lineRatioCatalogA->setLineAmplitude("CaII_t1A_8500.35_A", 0.00178621);
    lineRatioCatalogA->setLineAmplitude("CaII_t2A_8544.42_A", 0.00196559);
    lineRatioCatalogA->setLineAmplitude("CaII_t3A_8664.5_A", 0.00225234);
    lineRatioCatalogA->setLineAmplitude("P11A_8865.22_A", 0.00228397);
    lineRatioCatalogA->setLineAmplitude("P11_8865.22_E", 3.62e-18);
    lineRatioCatalogA->setLineAmplitude("P10A_9017.38_A", 0.00229739);
    lineRatioCatalogA->setLineAmplitude("P10_9017.38_E", 4.82e-18);
    lineRatioCatalogA->setLineAmplitude("P9_9231.55_E", 6.52e-18);
    lineRatioCatalogA->setLineAmplitude("P9A_9231.55_A", 0.00244926);
    lineRatioCatalogA->setLineAmplitude("P8A_9548.59_A", 0.00251125);
    lineRatioCatalogA->setLineAmplitude("P8_9548.59_E", 9.02e-18);
    lineRatioCatalogA->setLineAmplitude("P7_10052.13_E", 1.3e-17);
    lineRatioCatalogA->setLineAmplitude("P7A_10052.13_A", 0.0);
    lineRatioCatalogA->setLineAmplitude("P6_10941.09_E", 1.44e-17);
    lineRatioCatalogA->setLineAmplitude("P6A_10941.09_A", 0.00247716);
    lineRatioCatalogA->setLineAmplitude("P5_12821.59_E", 3.2e-17);
    lineRatioCatalogA->setLineAmplitude("P5A_12821.59_A", 0.00254503);
    lineRatioCatalogA->addVelocity("em_vel", 420.0);
    lineRatioCatalogA->addVelocity("abs_vel", 700.0);
    lineRatioCatalogA->setAsymProfileAndParams("ASYMFIXED",
                                               TAsymParams(2.0, 1.5, 0.0));
    lineRatioCatalogA->setPrior(0.5);
    lineRatioCatalogA->setIsmIndex(0);
  }
  std::shared_ptr<CLineRatioCatalog> lineRatioCatalogA =
      std::make_shared<CLineRatioCatalog>(
          "ratioCatalogA", *fixture_FullLineCatalog().lineCatalog);
};
 */
/* class fixture_LineRatioCatalogB {
public:
  fixture_LineRatioCatalogB() {

    lineRatioCatalogB->setLineAmplitude("LyGE_972.53_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("LyGA_972.53_A", 0.747995);
    lineRatioCatalogB->setLineAmplitude("LyBA_1025.72_A", 0.545877);
    lineRatioCatalogB->setLineAmplitude("LyBE_1025.72_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("OVI1031E_1031.93_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("OVI1031A_1031.93_A", 0.285283);
    lineRatioCatalogB->setLineAmplitude("OVI1037E_1037.62_E", 9.58557e-16);
    lineRatioCatalogB->setLineAmplitude("OVI1037A_1037.62_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("CIII1176a_1175.71_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("CIII1176b_1176.37_A", 0.254281);
    lineRatioCatalogB->setLineAmplitude("SiII1190_1190.42_A", 0.238655);
    lineRatioCatalogB->setLineAmplitude("SiII1193_1193.29_A", 0.103038);
    lineRatioCatalogB->setLineAmplitude("SiIII1206_1206.5_A", 0.227819);
    lineRatioCatalogB->setLineAmplitude("LyAA_1216.03_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("LyAE_1216.03_E", 1.26039e-14);
    lineRatioCatalogB->setLineAmplitude("SiII1260_1260.42_A", 0.206711);
    lineRatioCatalogB->setLineAmplitude("OI1302_1302.15_A", 0.213471);
    lineRatioCatalogB->setLineAmplitude("CII_1334.53_A", 0.147257);
    lineRatioCatalogB->setLineAmplitude("NiII_1370.5_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("SiIV1393_1393.74_A", 0.169026);
    lineRatioCatalogB->setLineAmplitude("SiIV1402_1402.77_A", 0.0390234);
    lineRatioCatalogB->setLineAmplitude("SiII1526_1526.71_A", 0.0970955);
    lineRatioCatalogB->setLineAmplitude("CIV1548A_1548.19_A", 0.0587184);
    lineRatioCatalogB->setLineAmplitude("CIV1548_1548.61_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("CIV1550A_1550.79_A", 0.0528466);
    lineRatioCatalogB->setLineAmplitude("CIV1550_1551.19_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("CI_1560.73_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("FeII1608_1608.44_A", 0.062399);
    lineRatioCatalogB->setLineAmplitude("HeII(doublet-2)_1640.33_E",
                                        2.76741e-16);
    lineRatioCatalogB->setLineAmplitude("HeII(doublet-3)_1640.47_E",
                                        4.15112e-16);
    lineRatioCatalogB->setLineAmplitude("OIII1661_1661.26_E", 2.56016e-16);
    lineRatioCatalogB->setLineAmplitude("OIII1666_1666.6_E", 3.42502e-16);
    lineRatioCatalogB->setLineAmplitude("AlII1670_1670.77_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("AlIII1854_1855.22_A", 0.0672517);
    lineRatioCatalogB->setLineAmplitude("AlIII1862_1863.29_A", 0.0654782);
    lineRatioCatalogB->setLineAmplitude("[CIII]1907_1906.68_E", 7.95165e-16);
    lineRatioCatalogB->setLineAmplitude("[CIII]1909_1908.73_E", 7.95165e-16);
    lineRatioCatalogB->setLineAmplitude("FeII2249_2250.49_A", 0.0328772);
    lineRatioCatalogB->setLineAmplitude("FeII2260_2261.39_A", 0.0673244);
    lineRatioCatalogB->setLineAmplitude("FeII2344_2344.84_A", 0.251161);
    lineRatioCatalogB->setLineAmplitude("FeII2374_2375.1_A", 0.150923);
    lineRatioCatalogB->setLineAmplitude("FeII2382_2383.4_A", 0.2111);
    lineRatioCatalogB->setLineAmplitude("FeII2586_2587.35_A", 0.207678);
    lineRatioCatalogB->setLineAmplitude("FeII2600_2600.87_A", 0.254946);
    lineRatioCatalogB->setLineAmplitude("MgII2796_2797.11_A", 0.199203);
    lineRatioCatalogB->setLineAmplitude("MgII_2799.12_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("MgII2803_2804.29_A", 0.0993576);
    lineRatioCatalogB->setLineAmplitude("MgI2852_2853.73_A", 0.0465985);
    lineRatioCatalogB->setLineAmplitude("HeI3190_3191.78_E", 9.44768e-18);
    lineRatioCatalogB->setLineAmplitude("[NeVb]_3346.81_E",
                                        9.785010000000001e-18);
    lineRatioCatalogB->setLineAmplitude("[NeVa]_3426.73_E", 3.85403e-17);
    lineRatioCatalogB->setLineAmplitude("FeI_3582.19_A", 0.0152131);
    lineRatioCatalogB->setLineAmplitude("[OII]3726_3727.09_E", 1.1e-15);
    lineRatioCatalogB->setLineAmplitude("[OII]3729_3729.88_E", 1.1e-15);
    lineRatioCatalogB->setLineAmplitude("FeVII_3758.46_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("H11_3771.65_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("H11A_3771.65_A", 0.164514);
    lineRatioCatalogB->setLineAmplitude("FeII3785_3785.47_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("H10_3798.93_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("H10A_3798.93_A", 0.191737);
    lineRatioCatalogB->setLineAmplitude("H9A_3836.43_A", 0.241925);
    lineRatioCatalogB->setLineAmplitude("H9_3836.43_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("NeIII_3869.05_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("H8_3890.11_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("H8A_3890.11_A", 0.184786);
    lineRatioCatalogB->setLineAmplitude("CaII_K_3934.73_A", 0.138547);
    lineRatioCatalogB->setLineAmplitude("CaII_H_3969.55_A", 0.132523);
    lineRatioCatalogB->setLineAmplitude("Hepsilon_3971.15_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("HdeltaA_4102.81_A", 0.137788);
    lineRatioCatalogB->setLineAmplitude("Hdelta_4102.81_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("GBand_4304.57_A", 0.0790501);
    lineRatioCatalogB->setLineAmplitude("HgammaA_4341.58_A", 0.0267463);
    lineRatioCatalogB->setLineAmplitude("Hgamma_4341.58_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("HbetaA_4862.72_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("Hbeta_4862.72_E", 4.42701e-16);
    lineRatioCatalogB->setLineAmplitude("[OIII](doublet-1/3)_4960.29_E",
                                        1.66e-15);
    lineRatioCatalogB->setLineAmplitude("[OIII](doublet-1)_5008.24_E", 5e-15);
    lineRatioCatalogB->setLineAmplitude("MgI5175_5176.71_A", 0.0479005);
    lineRatioCatalogB->setLineAmplitude("HeI5876_5877.41_E", 5.05679e-18);
    lineRatioCatalogB->setLineAmplitude("NaD_5895.6_A", 0.0281315);
    lineRatioCatalogB->setLineAmplitude("[OI]6301_6303.05_E", 9.86107e-17);
    lineRatioCatalogB->setLineAmplitude("[NII](doublet-1/2.95)_6549.86_E",
                                        7.03351e-17);
    lineRatioCatalogB->setLineAmplitude("HalphaA_6564.61_A", 0.0);
    lineRatioCatalogB->setLineAmplitude("Halpha_6564.61_E", 1.89434e-15);
    lineRatioCatalogB->setLineAmplitude("[NII](doublet-1)_6585.27_E",
                                        2.07489e-16);
    lineRatioCatalogB->setLineAmplitude("[SII]6716_6718.29_E", 4.15341e-16);
    lineRatioCatalogB->setLineAmplitude("[SII]6731_6732.68_E", 3.19644e-16);
    lineRatioCatalogB->setLineAmplitude("ArIII7136_7138.73_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("OII7319_7321.0_E", 0.0);
    lineRatioCatalogB->setLineAmplitude("OII7330_7332.2_E", 4.25514e-17);
    lineRatioCatalogB->setLineAmplitude("ArIII7751_7753.2_E", 2.72805e-17);
    lineRatioCatalogB->setLineAmplitude("CaII_t1A_8500.35_A", 0.00260494);
    lineRatioCatalogB->setLineAmplitude("CaII_t2A_8544.42_A", 0.00284238);
    lineRatioCatalogB->setLineAmplitude("CaII_t3A_8664.5_A", 0.00296519);
    lineRatioCatalogB->setLineAmplitude("P11A_8865.22_A", 0.00299126);
    lineRatioCatalogB->setLineAmplitude("P11_8865.22_E", 1.71e-17);
    lineRatioCatalogB->setLineAmplitude("P10A_9017.38_A", 0.00301032);
    lineRatioCatalogB->setLineAmplitude("P10_9017.38_E", 2.28e-17);
    lineRatioCatalogB->setLineAmplitude("P9_9231.55_E", 3.09e-17);
    lineRatioCatalogB->setLineAmplitude("P9A_9231.55_A", 0.00303573);
    lineRatioCatalogB->setLineAmplitude("P8A_9548.59_A", 0.00307068);
    lineRatioCatalogB->setLineAmplitude("P8_9548.59_E", 4.27e-17);
    lineRatioCatalogB->setLineAmplitude("P7_10052.13_E", 6.17e-17);
    lineRatioCatalogB->setLineAmplitude("P7A_10052.13_A", 0.0031203);
    lineRatioCatalogB->setLineAmplitude("P6_10941.09_E", 6.65e-17);
    lineRatioCatalogB->setLineAmplitude("P6A_10941.09_A", 0.00319345);
    lineRatioCatalogB->setLineAmplitude("P5_12821.59_E", 1.52e-16);
    lineRatioCatalogB->setLineAmplitude("P5A_12821.59_A", 0.00330484);
    lineRatioCatalogB->addVelocity("em_vel", 500.0);
    lineRatioCatalogB->addVelocity("abs_vel", 620.0);
    lineRatioCatalogB->setAsymProfileAndParams("ASYMFIXED",
                                               TAsymParams(2.0, 1.5, 1.0));
    lineRatioCatalogB->setPrior(0.5);
  }
  std::shared_ptr<CLineRatioCatalog> lineRatioCatalogB =
      std::make_shared<CLineRatioCatalog>(
          "ratioCatalogB", *fixture_FullLineCatalog().lineCatalog);
};
 */
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

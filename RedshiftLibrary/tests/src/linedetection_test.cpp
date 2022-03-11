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
#include "RedshiftLibrary/operator/linedetection.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/gaussianfit/gaussianfit.h"
#include "RedshiftLibrary/operator/linedetection.h"
#include "RedshiftLibrary/operator/linedetectionresult.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/common/exception.h"

#include <limits>
#include <iostream>
#include <cmath>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(test_linedetection)

BOOST_AUTO_TEST_CASE(XMadFind){
  CLineDetection lineDetection = CLineDetection( CLine::nType_Emission,0.5,0.6,0.7,0.8,0.9, true);

  TFloat64List x(10);
  x[0] = 5.;
  x[1] = 2.5;
  x[2] = 1.5;
  x[3] = 4.5;
  x[4] = 3.;

  Float64 returnValue = lineDetection.XMadFind(x.data(), 5 , 3.);
  BOOST_CHECK_CLOSE( returnValue, 1.5, 1e-12);

  x.resize(10);
  x[0] = 1.;
  x[1] = 1.;
  x[2] = 2.;
  x[3] = 3.;
  x[4] = 5.;
  x[5] = 3.;
  x[6] = 2.;
  x[7] = 1.;
  x[8] = 1.;
  x[9] = 1.;  
  
  CMedian<Float64> medianProcessor;
  Float64 med = medianProcessor.Find( x );
  BOOST_CHECK_CLOSE( med, 1.5, 1e-12);

  Float64 xmed = lineDetection.XMadFind(x.data(), 5 , med);
  BOOST_CHECK_CLOSE( xmed, 0.5, 1e-12);
}

BOOST_AUTO_TEST_CASE(ComputeFluxes){
  CLineDetection lineDetection = CLineDetection( CLine::nType_Emission,0.5,0.6,0.7,0.8,0.9, true);
  CSpectrum spc = CSpectrum();
  CSpectrumSpectralAxis spectralAxis(10, false );
  Float64* spcAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<spectralAxis.GetSamplesCount(); k++){
    spcAxis[k]=k;
  }

  CSpectrumFluxAxis modelfluxAxis(10);
  Float64* modelSamples = modelfluxAxis.GetSamples();
  modelSamples[0] = 1.;
  modelSamples[1] = 1.;
  modelSamples[2] = 2.;
  modelSamples[3] = 3.;
  modelSamples[4] = 5.;
  modelSamples[5] = 3.;
  modelSamples[6] = 2.;
  modelSamples[7] = 1.;
  modelSamples[8] = 1.;
  modelSamples[9] = 1.;

  spc.SetSpectralAndFluxAxes(std::move(spectralAxis), modelfluxAxis);
  Float64 winsize = 10000.;
  TInt32Range range = TInt32Range(0,9);
  TFloat64List mask = TFloat64List();

  Float64 maxFluxnoContinuum = 5.;
  Float64 noise = 12.;

  Float64 ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise);

  BOOST_CHECK_CLOSE( ratioAmp, 7.0, 1e-12);
  BOOST_CHECK_CLOSE( noise, 0.5, 1e-12);
  BOOST_CHECK_CLOSE( maxFluxnoContinuum, 3.5, 1e-12);

  winsize = 4.;
  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 2.0, 1e-12);
  BOOST_CHECK_CLOSE( noise, 1.0, 1e-12);
  BOOST_CHECK_CLOSE( maxFluxnoContinuum, 2.0, 1e-12);

  winsize = 10000.;
  mask.resize(10);
  for(Int32 k=0; k<10; k++){
    mask[k]=1.;
  }
  mask[4] = 0.;
  mask[5] = 0.;

  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise);
  BOOST_CHECK( ratioAmp == std::numeric_limits<double>::infinity() );  // return is inf because of noise = 0
  BOOST_CHECK_CLOSE( noise, 0.0, 1e-12);
  BOOST_CHECK_CLOSE( maxFluxnoContinuum, 2.0, 1e-12);

  modelSamples[8] = 1.2;
  modelSamples[9] = 1.1;
  spc.SetFluxAxis(modelfluxAxis);

  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 37./3, 1e-12);
  BOOST_CHECK_CLOSE( noise, 0.15, 1e-12);
  BOOST_CHECK_CLOSE( maxFluxnoContinuum, 1.85, 1e-12);

  mask[4] = 1.;
  mask[5] = 1.;
  range = TInt32Range(0,3);
  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 3.0, 1e-12);
  BOOST_CHECK_CLOSE( noise, 0.5, 1e-12);
  BOOST_CHECK_CLOSE( maxFluxnoContinuum, 1.5, 1e-12);

  range = TInt32Range(0,9);
  CSpectrumNoiseAxis& error = modelfluxAxis.GetError();
  error[0] = 0.5;
  error[1] = 0.3;
  error[2] = 0.8;
  error[3] = 0.8;
  error[4] = 0.9;
  error[5] = 0.7;
  error[6] = 0.6;
  error[7] = 0.8;
  error[8] = 0.9;
  error[9] = 0.3;

  spc.SetFluxAxis(modelfluxAxis);
  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 3.4/0.66, 1e-12);
  BOOST_CHECK_CLOSE( noise, 0.6, 1e-12);
  BOOST_CHECK_CLOSE( maxFluxnoContinuum, 3.4, 1e-12);

  range = TInt32Range(0,10);
  BOOST_CHECK_THROW(lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise), GlobalException);

  range = TInt32Range(-2,9);
  BOOST_CHECK_THROW(lineDetection.ComputeFluxes(spc, winsize, range, mask, &maxFluxnoContinuum, &noise), GlobalException);
}

BOOST_AUTO_TEST_CASE(RemoveStrongFromSpectra){
  CLineDetection lineDetection = CLineDetection( CLine::nType_Emission,0.5,0.6,0.7,0.8,0.9, true);
  Int32 n = 200;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* spcAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    spcAxis[k]=k;
  }
  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(n);
  for(Int32 k=0; k<n; k++){
    modelfluxAxis[k]=k;
  }

  Float64 sigma1 = 0.3;
  Float64 mu1 = 20.;
  Float64 A1 = 1.2;
  for(Int32 k=mu1-10; k<=mu1+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1));
  }
  CLineProfile_ptr profilesym{std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM()) };
  CLine line1 = CLine("Line1",mu1, 2, profilesym->Clone(), 2, A1, sigma1, 5.6);

  Float64 sigma2 = 0.5;
  Float64 mu2 = 70.;
  Float64 A2 = 2.2;
  for(Int32 k=mu2-10; k<=mu2+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1));
  }
  CLine line2 = CLine("Line2",mu2, 2, profilesym->Clone(), 2, A2, sigma2, 5.8);

  CSpectrum spc = CSpectrum(std::move(spectralAxis),std::move(modelfluxAxis));

  CLineDetectionResult lineDetectionResult;
  CLineCatalog::TLineVector strongLines;
  strongLines.push_back(line1);
  strongLines.push_back(line2);

  TInt32RangeList selectedretestPeaks;
  selectedretestPeaks.push_back(TInt32Range(5,35));
  selectedretestPeaks.push_back(TInt32Range(60,80));

  CLineDetection::TGaussParamsList selectedgaussparams;
  selectedgaussparams.push_back(CLineDetection::SGaussParams(2.,0.2,0.3));
  selectedgaussparams.push_back(CLineDetection::SGaussParams(7.,1.,0.5));
  Float64 winsize = 100.;
  Float64 cut = -3.;

  lineDetection.RemoveStrongFromSpectra(spc, lineDetectionResult, strongLines, selectedretestPeaks, selectedgaussparams, winsize, cut);
  BOOST_CHECK_EQUAL(lineDetectionResult.LineCatalog.GetList().size(), 2);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetCut(), 1.7505788080267244, 1e-12);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetCut(), 2, 1e-12);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetAmplitude(), 0.3, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetAmplitude(), 0.5, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetPosition(), 2.0, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetPosition(), 7.0, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetWidth(), 0.2, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetWidth(), 1.0, 1e-6);
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[0].GetProfile().GetName() == profilesym->GetName());
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[1].GetProfile().GetName() == profilesym->GetName());
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[0].GetIsStrong() == false);
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[1].GetIsStrong() == false);
}

BOOST_AUTO_TEST_CASE(Retest){
  CLineDetection lineDetection = CLineDetection( CLine::nType_Emission,0.5,0.6,0.7,0.8,0.9, true);
  Int32 n = 200;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* spcAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    spcAxis[k]=k;
  }
  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(n);
  for(Int32 k=0; k<n; k++){
    modelfluxAxis[k]=k;
  }
  Float64 sigma1 = 0.3;
  Float64 mu1 = 20.;
  Float64 A1 = 1.2;
  for(Int32 k=mu1-10; k<=mu1+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1));
  }
  CLineProfile_ptr profilesym{std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM()) };
  CLine line1 = CLine("Line1",mu1, 2, profilesym->Clone(), 2, A1, sigma1, 5.6);

  Float64 sigma2 = 0.5;
  Float64 mu2 = 70.;
  Float64 A2 = 2.2;
  for(Int32 k=mu2-10; k<=mu2+10; k++){
    modelfluxAxis[k]+=A2/(sigma2 *2.506597694086548) *exp(-(k-mu2)*(k-mu2)/(2*sigma2*sigma2));
  }
  CLine line2 = CLine("Line2",mu2, 2, profilesym->Clone(), 2, A2, sigma2, 5.8);

  CSpectrum spc = CSpectrum(std::move(spectralAxis),std::move(modelfluxAxis));

  CLineDetectionResult lineDetectionResult;
  CLineCatalog::TLineVector strongLines;
  strongLines.push_back(line1);
  strongLines.push_back(line2);

  TInt32RangeList retestPeaks;
  retestPeaks.push_back(TInt32Range(5,35));
  retestPeaks.push_back(TInt32Range(60,80));
  retestPeaks.push_back(TInt32Range(101,199));

  CLineDetection::TGaussParamsList selectedgaussparams;
  selectedgaussparams.push_back(CLineDetection::SGaussParams(2.,0.2,0.3));
  selectedgaussparams.push_back(CLineDetection::SGaussParams(7.,1.,0.5));
  Float64 winsize = 100.;
  Float64 cut = -3.;

  lineDetection.Retest(spc, lineDetectionResult, retestPeaks, selectedgaussparams, strongLines, winsize, cut);
  BOOST_CHECK_EQUAL(lineDetectionResult.LineCatalog.GetList().size() , 2);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetCut(), 1.7505788080267244, 1e-12);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetCut(), 1.6729987967017514, 1e-12);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetAmplitude(), 0.3, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetAmplitude(), 0.5, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetPosition(), 2.0, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetPosition(), 7.0, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetWidth(), 0.2, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[1].GetWidth(), 1.0, 1e-6);
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[0].GetProfile().GetName() == profilesym->GetName());
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[1].GetProfile().GetName() == profilesym->GetName());
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[0].GetIsStrong() == false);
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[1].GetIsStrong() == false);

  TInt32RangeList retestPeaks2;
  retestPeaks2.push_back(TInt32Range(5,35));

  modelfluxAxis = CSpectrumFluxAxis( spc.GetSampleCount());
  for(Int32 k=0; k<spc.GetSampleCount(); k++){
    modelfluxAxis[k]=k;
  }
  for(Int32 k=mu1-10; k<=mu1+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1));
  }
  for(Int32 k=mu2-10; k<=mu2+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1));
  }
  spc.SetFluxAxis(std::move(modelfluxAxis));

  CLineDetectionResult lineDetectionResult2;
  lineDetection.Retest(spc, lineDetectionResult2, retestPeaks2, selectedgaussparams, strongLines, winsize, cut);
  BOOST_CHECK_EQUAL(lineDetectionResult2.LineCatalog.GetList().size(), 1);
  BOOST_CHECK_CLOSE(lineDetectionResult2.LineCatalog.GetList()[0].GetCut(), 1.7505788080267244, 1e-12);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetAmplitude(), 0.3, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetPosition(), 2.0, 1e-6);
  BOOST_CHECK_CLOSE(lineDetectionResult.LineCatalog.GetList()[0].GetWidth(), 0.2, 1e-6);
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[0].GetProfile().GetName() == profilesym->GetName());
  BOOST_CHECK(lineDetectionResult.LineCatalog.GetList()[0].GetIsStrong() == false);
}

BOOST_AUTO_TEST_CASE(LimitGaussianFitStartAndStop){
  Int32 n = 250;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* fluxAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    fluxAxis[k]=k;
  }

  TInt32RangeList peak;
  peak.push_back(TInt32Range(-5,35));
  peak.push_back(TInt32Range(30,70));
  peak.push_back(TInt32Range(80,120));
  peak.push_back(TInt32Range(160,200));
  peak.push_back(TInt32Range(190,230));

  CLineDetection lineDetection = CLineDetection( CLine::nType_Emission,0.5,0.6,0.7,0.8,40, true);
  TInt32Range range = lineDetection.LimitGaussianFitStartAndStop( 0,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 0);
  BOOST_CHECK_EQUAL(range.GetEnd(), 30);

  range = lineDetection.LimitGaussianFitStartAndStop( 1,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 35);
  BOOST_CHECK_EQUAL(range.GetEnd(), 70);

  range = lineDetection.LimitGaussianFitStartAndStop( 2,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 80);
  BOOST_CHECK_EQUAL(range.GetEnd(), 120);

  range = lineDetection.LimitGaussianFitStartAndStop( 3,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 160);
  BOOST_CHECK_EQUAL(range.GetEnd(), 190);

  range = lineDetection.LimitGaussianFitStartAndStop( 4,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 200);
  BOOST_CHECK_EQUAL(range.GetEnd(), 230);
}

void addLine(CSpectrumFluxAxis& spectrumFluxAxis, Float64 sigma, Float64 mu, Float64 A){
  for(Int32 k=mu-sigma*5; k<=mu+sigma*5; k++){
    spectrumFluxAxis[k]+= A*exp(-(k-mu)*(k-mu)/(2*sigma*sigma));
  }
}

BOOST_AUTO_TEST_CASE(Compute){
  CLineDetection lineDetection = CLineDetection(CLine::nType_Emission);

  Int32 n = 2000;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* spcAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    spcAxis[k]=k;
  }
  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(n);
  for(Int32 k=0; k<n; k++){
    modelfluxAxis[k]=k*0.0001;
  }
  CSpectrumNoiseAxis& error = modelfluxAxis.GetError();
  for(Int32 k=0; k<n; k++){
    error[k]=k*0.0001;
  }
  TInt32RangeList resPeaks;

  addLine(modelfluxAxis,4.,40.,1.5);
  resPeaks.push_back(TInt32Range(20,60));

  addLine(modelfluxAxis,0.5,100.,1.5);
  resPeaks.push_back(TInt32Range(60,140));

  addLine(modelfluxAxis,70.,500.,1.5);
  resPeaks.push_back(TInt32Range(150,850));

  addLine(modelfluxAxis,4.,1000.,1.5);
  addLine(modelfluxAxis,4.,1030.,1.5);
  resPeaks.push_back(TInt32Range(950,1100));

  addLine(modelfluxAxis,40.,1200.,3.5);
  addLine(modelfluxAxis,4.,1140.,1.5);
  resPeaks.push_back(TInt32Range(1130,1150));

  addLine(modelfluxAxis,4.,1450.,-13.5);
  resPeaks.push_back(TInt32Range(1400,1500));
  resPeaks.push_back(TInt32Range(1500,1500));

  // for(Int32 k=1100 ;k<=1500; k++){
  //   BOOST_TEST_MESSAGE("modelfluxAxis "<< modelfluxAxis[k]);
  // }
  CSpectrum spc = CSpectrum(std::move(spectralAxis),std::move(modelfluxAxis));


  TLambdaRange lambdaRange = TLambdaRange(0.0,2000.0); //useless


  std::shared_ptr<const CLineDetectionResult> res = lineDetection.Compute(spc, lambdaRange, resPeaks, resPeaks);

  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};

  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[0], "Peak_0 : line detected successfully");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[1], "Peak_1 : fwhm<m_minsize");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[2], "Peak_2 : fwhm>m_maxsize");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[3], "Peak_3 : gaussAmp far from spectrum max_value");
  //bug here two status for one peak
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[4], "Peak_3 : gaussAmp far from spectrum max_value (Angstrom)");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[5], "Peak_4 : ratioAmp<m_cut (1.526190<5.000000)");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[6], "Peak_5 : GaussAmp negative");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[7], "Peak_6 : Fitting failed");

  BOOST_CHECK_EQUAL(res->LineCatalog.GetList().size(), 1);
  BOOST_CHECK_CLOSE(res->LineCatalog.GetList()[0].GetAmplitude(), 1.5, 1e-6);
  BOOST_CHECK_CLOSE(res->LineCatalog.GetList()[0].GetPosition(), 40.0, 1e-6);
  BOOST_CHECK_CLOSE(res->LineCatalog.GetList()[0].GetWidth(), 4.0, 1e-6);
  BOOST_CHECK(res->LineCatalog.GetList()[0].GetProfile().GetName() == profilesym->GetName());
  BOOST_CHECK(res->LineCatalog.GetList()[0].GetIsStrong() == true);
  BOOST_CHECK(res->LineCatalog.GetList()[0].GetIsEmission() == true);
}


BOOST_AUTO_TEST_SUITE_END()

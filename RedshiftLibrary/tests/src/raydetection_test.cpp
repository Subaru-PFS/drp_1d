#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/common/median.h>
#include <RedshiftLibrary/gaussianfit/gaussianfit.h>
#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>



#include <time.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>


using namespace NSEpic;
BOOST_AUTO_TEST_SUITE(test_raydetection)


//http://bloglitb.blogspot.fr/2010/07/access-to-private-members-thats-easy.html
// allow to acces to private method
template<typename Tag>
struct result {
  /* export it ... */
  typedef typename Tag::type type;
  static type ptr;
};

template<typename Tag>
typename result<Tag>::type result<Tag>::ptr;

template<typename Tag, typename Tag::type p>
struct rob : result<Tag> {
  /* fill it ... */
  struct filler {
    filler() {
      result<Tag>::ptr = p; }
  };
  static filler filler_obj;
};

template<typename Tag, typename Tag::type p>
typename rob<Tag, p>::filler rob<Tag, p>::filler_obj;


struct CLineDetectionXMadFind { typedef Float64(CLineDetection::*type)(const Float64* x, Int32 n, Float64 median ); };
template class rob<CLineDetectionXMadFind, &CLineDetection::XMadFind>;


struct CLineDetectionRemoveStrongFromSpectra { typedef bool(CLineDetection::*type)(const CSpectrum& spectrum, CLineDetectionResult& result,  CRayCatalog::TRayVector strongLines, TInt32RangeList selectedretestPeaks, CLineDetection::TGaussParamsList selectedgaussparams, Float64 winsize, Float64 cut); };
template class rob<CLineDetectionRemoveStrongFromSpectra, &CLineDetection::RemoveStrongFromSpectra>;


struct CLineDetectionRetest { typedef bool(CLineDetection::*type)( const CSpectrum& spectrum, CLineDetectionResult& result, TInt32RangeList retestPeaks,  CLineDetection::TGaussParamsList retestGaussParams, CRayCatalog::TRayVector strongLines, Int32 winsize, Float64 cut); };
template class rob<CLineDetectionRetest, &CLineDetection::Retest>;

struct CLineDetectionLimitGaussianFitStartAndStop { typedef TInt32Range(CLineDetection::*type)(Int32 i, const TInt32RangeList& peaksBorders, Int32 len, const CSpectrumSpectralAxis spectralAxis ); };
template class rob<CLineDetectionLimitGaussianFitStartAndStop, &CLineDetection::LimitGaussianFitStartAndStop>;

BOOST_AUTO_TEST_CASE(XMadFind){
  CLineDetection lineDetection = CLineDetection( 1,0.5,0.6,0.7,0.8,0.9, true);

  Float64* x = (Float64*) calloc( 5, sizeof( Float64 ) );
  x[0] = 5.;
  x[1] = 2.5;
  x[2] = 1.5;
  x[3] = 4.5;
  x[4] = 3.;

  Float64 returnValue = (lineDetection.*result<CLineDetectionXMadFind>::ptr)(x, 5 , 3.);
  BOOST_CHECK_CLOSE( returnValue, 1.5, 1e-12);



  x = (Float64*) calloc( 10, sizeof( Float64 ) );
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
  Float64 med = medianProcessor.Find( x, 10 );
  BOOST_CHECK_CLOSE( med, 1.5, 1e-12);
  Float64 xmed = (lineDetection.*result<CLineDetectionXMadFind>::ptr)(x, 5 , med);
  BOOST_CHECK_CLOSE( xmed, 0.5, 1e-12);
  free(x);
}

BOOST_AUTO_TEST_CASE(ComputeFluxes){
  CLineDetection lineDetection = CLineDetection( 1,0.5,0.6,0.7,0.8,0.9, true);
  CSpectrum spc =  CSpectrum();
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(10, false );
  Float64* fluxAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<spectralAxis.GetSamplesCount(); k++){
    fluxAxis[k]=k;
  }

  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(10);
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
  Float64 winsize = 10000;
  TInt32Range range = TInt32Range(0,10);
  TFloat64List mask= TFloat64List();

  Float64 maxFluxnoContinuum = 5.;
  Float64 noise = 12.;

  Float64 ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask,&maxFluxnoContinuum, &noise);

  BOOST_CHECK_CLOSE( ratioAmp, 7., 1e-12);
  BOOST_CHECK_CLOSE( noise , 0.5,1e-12);
  BOOST_CHECK_CLOSE(maxFluxnoContinuum, 3.5 , 1e-12);

  winsize = 4;
  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask,&maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 2., 1e-12);
  BOOST_CHECK_CLOSE( noise , 1,1e-12);
  BOOST_CHECK_CLOSE(maxFluxnoContinuum, 2 , 1e-12);

  winsize = 10000;
  mask.resize(10);
  for(Int32 k=0; k<10; k++){
    mask[k]=1;
  }
  mask[4] = 0;
  mask[5] = 0;

  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask,&maxFluxnoContinuum, &noise);
  //BOOST_CHECK_CLOSE( ratioAmp, 1, 1e-12);// return is inf because of noise = 0
  BOOST_CHECK_CLOSE( noise , 0,1e-12);
  BOOST_CHECK_CLOSE(maxFluxnoContinuum, 2 , 1e-12);

  modelSamples[8] = 1.2;
  modelSamples[9] = 1.1;
  spc.SetFluxAxis(modelfluxAxis);

  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask,&maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 37./3, 1e-12);
  BOOST_CHECK_CLOSE( noise , 0.15,1e-12);
  BOOST_CHECK_CLOSE(maxFluxnoContinuum, 1.85 , 1e-12);


  mask[4] = 1;
  mask[5] = 1;
  range = TInt32Range(0,3);
  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask,&maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 3.0, 1e-12);
  BOOST_CHECK_CLOSE( noise , 0.5,1e-12);
  BOOST_CHECK_CLOSE(maxFluxnoContinuum, 1.5 , 1e-12);


  range = TInt32Range(0,10);
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
  ratioAmp = lineDetection.ComputeFluxes(spc, winsize, range, mask,&maxFluxnoContinuum, &noise);
  BOOST_CHECK_CLOSE( ratioAmp, 3.4/0.66, 1e-12);
  BOOST_CHECK_CLOSE( noise , 0.6, 1e-12);
  BOOST_CHECK_CLOSE( maxFluxnoContinuum, 3.4 , 1e-12);

}
BOOST_AUTO_TEST_CASE(RemoveStrongFromSpectra){
  CLineDetection lineDetection = CLineDetection( 1,0.5,0.6,0.7,0.8,0.9, true);
  Int32 n = 200;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* fluxAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    fluxAxis[k]=k;
  }
  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(n);
  for(Int32 k=0; k<n; k++){
    modelfluxAxis[k]=k;
  }

  Float64 sigma1 = 0.3;
  Float64 mu1 = 20.;
  Float64 A1 = 1.2;
  //for(Int32 k=mu1-10; k<=mu1+10; k++){
  //  modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1)) ;
  for(Int32 k=mu1-10; k<mu1+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1)/(2*sigma1)) ;
  }
  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
  CRay ray1 = CRay("Ray1",mu1, 2, profilesym, 2, A1, sigma1, 5.6);

  Float64 sigma2 = 0.5;
  Float64 mu2 = 70.;
  Float64 A2 = 2.2;
  //for(Int32 k=mu2-10; k<=mu2+10; k++){
  //  modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1)) ;
  for(Int32 k=mu2-10; k<mu2+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1)/(2*sigma1)) ;
  }
  CRay ray2 = CRay("Ray2",mu2, 2, profilesym, 2, A2, sigma2, 5.8);

  CSpectrum spc = CSpectrum(std::move(spectralAxis),std::move(modelfluxAxis));

  CLineDetectionResult lineDetectionResult;
  CRayCatalog::TRayVector strongLines;
  strongLines.push_back(ray1);
  strongLines.push_back(ray2);

  TInt32RangeList selectedretestPeaks;
  selectedretestPeaks.push_back(TInt32Range(5,35));
  selectedretestPeaks.push_back(TInt32Range(60,80));

  CLineDetection::TGaussParamsList selectedgaussparams;
  selectedgaussparams.push_back(CLineDetection::SGaussParams(2.,0.2,0.3));
  selectedgaussparams.push_back(CLineDetection::SGaussParams(7.,1.,0.5));
  Float64 winsize =100;
  Float64 cut = -3;


  (lineDetection.*result<CLineDetectionRemoveStrongFromSpectra>::ptr)(spc, lineDetectionResult, strongLines, selectedretestPeaks, selectedgaussparams, winsize, cut);
  BOOST_CHECK_EQUAL(lineDetectionResult.RayCatalog.GetList().size() , 2);
  BOOST_CHECK_CLOSE(lineDetectionResult.RayCatalog.GetList()[0].GetCut() , 1.7594187467864928, 1e-12);
  BOOST_CHECK_CLOSE(lineDetectionResult.RayCatalog.GetList()[1].GetCut() , 2.0, 1e-12);

}

BOOST_AUTO_TEST_CASE(Retest){
  CLineDetection lineDetection = CLineDetection( 1,0.5,0.6,0.7,0.8,0.9, true);
  Int32 n = 200;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* fluxAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    fluxAxis[k]=k;
  }
  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(n);
  for(Int32 k=0; k<n; k++){
    modelfluxAxis[k]=k;
  }
  Float64 sigma1 = 0.3;
  Float64 mu1 = 20.;
  Float64 A1 = 1.2;
  //for(Int32 k=mu1-10; k<=mu1+10; k++){
  //  modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1)) ;
  for(Int32 k=mu1-10; k<mu1+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1)/(2*sigma1)) ;
  }
  std::shared_ptr<CLineProfile> profilesym{std::make_shared<CLineProfileSYM>()};
  CRay ray1 = CRay("Ray1",mu1, 2, profilesym, 2, A1, sigma1, 5.6);

  Float64 sigma2 = 0.5;
  Float64 mu2 = 70.;
  Float64 A2 = 2.2;
  //for(Int32 k=mu2-10; k<=mu2+10; k++){
  //  modelfluxAxis[k]+=A2/(sigma2 *2.506597694086548) *exp(-(k-mu2)*(k-mu2)/(2*sigma2*sigma2)) ;
  for(Int32 k=mu2-10; k<mu2+10; k++){
    modelfluxAxis[k]+=A2/(sigma2 *2.506597694086548) *exp(-(k-mu2)*(k-mu2)/(2*sigma2)/(2*sigma2)) ;
  }
  CRay ray2 = CRay("Ray2",mu2, 2, profilesym, 2, A2, sigma2, 5.8);

  CSpectrum spc = CSpectrum(std::move(spectralAxis),std::move(modelfluxAxis));


  CLineDetectionResult lineDetectionResult;
  CRayCatalog::TRayVector strongLines;
  strongLines.push_back(ray1);
  strongLines.push_back(ray2);

  TInt32RangeList retestPeaks;
  retestPeaks.push_back(TInt32Range(5,35));
  retestPeaks.push_back(TInt32Range(60,80));
  retestPeaks.push_back(TInt32Range(100,199));

  CLineDetection::TGaussParamsList selectedgaussparams;
  selectedgaussparams.push_back(CLineDetection::SGaussParams(2.,0.2,0.3));
  selectedgaussparams.push_back(CLineDetection::SGaussParams(7.,1.,0.5));
  Int32 winsize = 100;
  Float64 cut = -3;

  (lineDetection.*result<CLineDetectionRetest>::ptr)(spc, lineDetectionResult,  retestPeaks, selectedgaussparams,strongLines, winsize, cut);
  BOOST_CHECK_EQUAL(lineDetectionResult.RayCatalog.GetList().size() , 2);
  BOOST_CHECK_CLOSE(lineDetectionResult.RayCatalog.GetList()[0].GetCut(), 1.7594187467864928, 1e-12);
  BOOST_CHECK_CLOSE(lineDetectionResult.RayCatalog.GetList()[1].GetCut(), 1.5603039861799142, 1e-12);


  TInt32RangeList retestPeaks2;
  retestPeaks2.push_back(TInt32Range(5,35));


  modelfluxAxis = CSpectrumFluxAxis( spc.GetSampleCount());
  for(Int32 k=0; k<spc.GetSampleCount(); k++){
    modelfluxAxis[k]=k;
  }
  //for(Int32 k=mu1-10; k<=mu1+10; k++){
  //  modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1)) ;
  for(Int32 k=mu1-10; k<mu1+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1)/(2*sigma1)) ;
  }
  //for(Int32 k=mu2-10; k<=mu2+10; k++){
  //  modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1*sigma1)) ;
  for(Int32 k=mu2-10; k<mu2+10; k++){
    modelfluxAxis[k]+=A1/(sigma1 *2.506597694086548) *exp(-(k-mu1)*(k-mu1)/(2*sigma1)/(2*sigma1)) ;
  }
  spc.SetFluxAxis(std::move(modelfluxAxis));


  CLineDetectionResult lineDetectionResult2;
  (lineDetection.*result<CLineDetectionRetest>::ptr)(spc, lineDetectionResult2,  retestPeaks2, selectedgaussparams, strongLines,winsize, cut);
  BOOST_CHECK_EQUAL(lineDetectionResult2.RayCatalog.GetList().size() , 1);
  BOOST_CHECK_CLOSE(lineDetectionResult2.RayCatalog.GetList()[0].GetCut(), 1.7594187467864928, 1e-12);
}


BOOST_AUTO_TEST_CASE(LimitGaussianFitStartAndStop){
  Int32 n = 200;
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
  //peak.push_back(TInt32Range(160,250));==> bug

  CLineDetection lineDetection = CLineDetection( 1,0.5,0.6,0.7,0.8,40, true);
  TInt32Range range = (lineDetection.*result<CLineDetectionLimitGaussianFitStartAndStop>::ptr)( 0,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 0);
  BOOST_CHECK_EQUAL(range.GetEnd(), 30);

  range = (lineDetection.*result<CLineDetectionLimitGaussianFitStartAndStop>::ptr)( 1,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 35);
  BOOST_CHECK_EQUAL(range.GetEnd(), 70);

  range = (lineDetection.*result<CLineDetectionLimitGaussianFitStartAndStop>::ptr)( 2,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 80);
  BOOST_CHECK_EQUAL(range.GetEnd(), 120);

  range = (lineDetection.*result<CLineDetectionLimitGaussianFitStartAndStop>::ptr)( 3,  peak, n, spectralAxis);
  BOOST_CHECK_EQUAL(range.GetBegin(), 160);
  BOOST_CHECK_EQUAL(range.GetEnd(), 199);

}

void addRay(CSpectrumFluxAxis& spectrumFluxAxis, Float64 sigma, Float64 mu, Float64 A){
  for(Int32 k=mu-sigma*5; k<=mu+sigma*5; k++){
    spectrumFluxAxis[k]+= A*exp(-(k-mu)*(k-mu)/(2*sigma*sigma));
  }
}

BOOST_AUTO_TEST_CASE(Compute){
  CLineDetection lineDetection = CLineDetection(CRay::nType_Emission);

  Int32 n = 2000;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* fluxAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    fluxAxis[k]=k;
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

  addRay(modelfluxAxis,4.,40.,1.5);
  resPeaks.push_back(TInt32Range(20,60));

  addRay(modelfluxAxis,0.5,100.,1.5);
  resPeaks.push_back(TInt32Range(60,140));

  addRay(modelfluxAxis,70.,500.,1.5);
  resPeaks.push_back(TInt32Range(150,850));

  addRay(modelfluxAxis,4.,1000.,1.5);
  addRay(modelfluxAxis,4.,1030.,1.5);
  resPeaks.push_back(TInt32Range(950,1100));

  addRay(modelfluxAxis,40.,1200.,3.5);
  addRay(modelfluxAxis,4.,1140.,1.5);
  resPeaks.push_back(TInt32Range(1130,1150));

  addRay(modelfluxAxis,4.,1450.,-13.5);
  resPeaks.push_back(TInt32Range(1400,1500));
  resPeaks.push_back(TInt32Range(1500,1500));

  // for(Int32 k=1100 ;k<=1500; k++){
  //   BOOST_TEST_MESSAGE("modelfluxAxis "<< modelfluxAxis[k]);
  // }
  CSpectrum spc = CSpectrum(std::move(spectralAxis),std::move(modelfluxAxis));


  TLambdaRange lambdaRange = TLambdaRange(0.0,2000.0); //useless


  std::shared_ptr<const CLineDetectionResult> res = lineDetection.Compute(spc, lambdaRange, resPeaks, resPeaks);


  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[0], "Peak_0 : line detected successfully");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[1], "Peak_1 : fwhm<m_minsize");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[2], "Peak_2 : fwhm>m_maxsize");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[3], "Peak_3 : gaussAmp far from spectrum max_value");
  //bug here two status for one peak
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[4], "Peak_3 : gaussAmp far from spectrum max_value (Angstrom)");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[5], "Peak_4 : ratioAmp<m_cut (1.526190<5.000000)");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[6], "Peak_5 : GaussAmp negative");
  BOOST_CHECK_EQUAL(res->PeakListDetectionStatus[7], "Peak_6 : Fitting failed");

  // bogus tests
  /*
  BOOST_CHECK_EQUAL(res->RayCatalog.GetList().size(), 1);
  BOOST_CHECK_CLOSE(res->RayCatalog.GetList()[0].GetAmplitude(), 1.5, 1e-5);
  BOOST_CHECK_CLOSE(res->RayCatalog.GetList()[0].GetPosition(),40, 1e-5);
  */
}


BOOST_AUTO_TEST_SUITE_END()

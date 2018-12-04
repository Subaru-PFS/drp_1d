#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/common/median.h>
#include <RedshiftLibrary/gaussianfit/gaussianfit.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>



#include <time.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <boost/test/unit_test.hpp>


using namespace NSEpic;
BOOST_AUTO_TEST_SUITE(gaussianfit_test)


void addRay(CSpectrumFluxAxis& spectrumFluxAxis , Float64 sigma, Float64 mu, Float64 A){
  for(Int32 k= mu-sigma*5; k<mu+sigma*5; k++){
    spectrumFluxAxis[k]+=A *exp(-(k-mu)*(k-mu)/(2*sigma)/(2*sigma)) ;
  }
}

Float64 precision = 1e-4;

BOOST_AUTO_TEST_CASE(GaussianFit){
  CSpectrum spc =  CSpectrum();

  Int32 n = 1000;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  Float64* fluxAxis = spectralAxis.GetSamples();
  for(Int32 k=0; k<n; k++){
    fluxAxis[k]=k;
  }
  spc.GetSpectralAxis() = spectralAxis;
  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(n);
  for(Int32 k=0; k<n; k++){
    modelfluxAxis[k]=k*0.000;
  }
  TFloat64List& error = modelfluxAxis.GetError();
  for(Int32 k=0; k<n; k++){
    error[k]=0.5;
  }
  addRay(modelfluxAxis, 4,40,1.5);
  // addRay(modelfluxAxis, 4,80,-1.5);

  modelfluxAxis[80]= 1.5;
  modelfluxAxis[81]= 0.1;
  modelfluxAxis[82]= 8.2;
  modelfluxAxis[83]= 3.5;
  modelfluxAxis[84]= 7.4;
  modelfluxAxis[85]= 6.2;
  modelfluxAxis[86]= 8.1;
  modelfluxAxis[87]= 6.3;
  modelfluxAxis[88]= 0.8;
  modelfluxAxis[89]= 1.4;
  modelfluxAxis[90]= 7.4;
  modelfluxAxis[91]= 5.2;
  modelfluxAxis[92]= 9.2;
  modelfluxAxis[93]= 4.2;
  modelfluxAxis[94]= 2.5;
  modelfluxAxis[95]= 3.3;
  modelfluxAxis[96]= 3.9;
  modelfluxAxis[97]= 8.7;
  modelfluxAxis[98]= 5.4;
  modelfluxAxis[99]= 1.2;


  addRay(modelfluxAxis, 4,150.5,1.5);
  spc.GetFluxAxis() = modelfluxAxis;

  CGaussianFit fitter;

  Float64 gaussAmp;
  Float64 gaussPos;
  Float64 gaussWidth;

  Float64 gaussAmpErr;
  Float64 gaussPosErr;
  Float64 gaussWidthErr;
  Float64 coeff0;
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  BOOST_CHECK_CLOSE(gaussAmp, 0, precision);
  BOOST_CHECK_CLOSE(gaussPos, 0., precision);
  BOOST_CHECK_CLOSE(gaussWidth, 0, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 0., precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 0, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 0., precision);

  CGaussianFit::EStatus status = fitter.Compute( spc, TInt32Range( 0, 60 ) );
  BOOST_CHECK_EQUAL(status, CGaussianFit::EStatus::nStatus_Success);
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  fitter.GetResultsPolyCoeff0( coeff0 );
  BOOST_CHECK_CLOSE(gaussAmp, 1.5001791032853808, precision);
  BOOST_CHECK_CLOSE(gaussPos, 39.999667102453245, precision);
  BOOST_CHECK_CLOSE(gaussWidth, 5.6577322718280509, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 0.48499052812117621, precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 1.8868145649073094, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 2.3664087587709757, precision);
  BOOST_CHECK_CLOSE(coeff0, -0.00021929044444247408, precision);


  status = fitter.Compute( spc, TInt32Range( 50, 50 ) );
  BOOST_CHECK_EQUAL(status,CGaussianFit::EStatus::nStatus_IllegalInput);

  status = fitter.Compute( spc, TInt32Range( 50, 52 ) );
  BOOST_CHECK_EQUAL(status,CGaussianFit::EStatus::nStatus_IllegalInput);


  status = fitter.Compute( spc, TInt32Range( 80, 100 ) );
  BOOST_CHECK_EQUAL(status, CGaussianFit::EStatus::nStatus_Success);
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  fitter.GetResultsPolyCoeff0( coeff0 );
  BOOST_CHECK_CLOSE(gaussAmp, -6151.3675460026689, precision);
  BOOST_CHECK_CLOSE(gaussPos, 90.13718519126094, precision);
  BOOST_CHECK_CLOSE(gaussWidth, 27.549647405642968, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 1295407.1033151825, precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 1.051122999448332, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 1472.7234076202869, precision);
  BOOST_CHECK_CLOSE(coeff0, 6155.824288046083, precision);


  status = fitter.Compute( spc, TInt32Range( 140, 160 ) );
  BOOST_CHECK_EQUAL(status, CGaussianFit::EStatus::nStatus_Success);
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  fitter.GetResultsPolyCoeff0( coeff0 );
  BOOST_CHECK_CLOSE(gaussAmp, 1.5000000000000018, precision);
  BOOST_CHECK_CLOSE(gaussPos, 150.5, precision);
  BOOST_CHECK_CLOSE(gaussWidth, 5.6568542494923824, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 38.966662308077282, precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 5.4530518026175931, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 53.959488049892876, precision);
  BOOST_CHECK_CLOSE(coeff0, -1.7933884369275193e-15, precision);


  // status = fitter.Compute( spc, TInt32Range( 60, 100 ) );
  // fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  // BOOST_CHECK_CLOSE(gaussAmp, 0.14963071532336933, 1e-6);
  // BOOST_CHECK_CLOSE(gaussPos, 80, 1e-6);

}

BOOST_AUTO_TEST_SUITE_END()

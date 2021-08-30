#include "RedshiftLibrary/operator/raydetection.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/gaussianfit/gaussianfit.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"

#include <iostream>
#include <cmath>
#include <boost/test/unit_test.hpp>


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(gaussianfit_test)

void addRay(CSpectrumFluxAxis& spectrumFluxAxis, Float64 sigma, Float64 mu, Float64 A){
  for(Int32 k=mu-sigma*5; k<=mu+sigma*5; k++){
    spectrumFluxAxis[k]+=A*exp(-(k-mu)*(k-mu)/(2*sigma*sigma));
  }
}

Float64 precision = 1e-12;

BOOST_AUTO_TEST_CASE(GaussianFit){

  Int32 n = 1000;
  CSpectrumSpectralAxis spectralAxis = CSpectrumSpectralAxis(n, false );
  for(Int32 k=0; k<n; k++){
    spectralAxis[k]=k;
  }

  CSpectrumFluxAxis modelfluxAxis = CSpectrumFluxAxis(n);
  for(Int32 k=0; k<n; k++){
    modelfluxAxis[k]=0;
  }
  CSpectrumNoiseAxis& error = modelfluxAxis.GetError();
  for(Int32 k=0; k<n; k++){
    error[k]=0.5;
  }

  addRay(modelfluxAxis, 4.,40.,1.5);
  addRay(modelfluxAxis, 4.,80.,4.5);
  addRay(modelfluxAxis, 3.2,120.85,-5.0);
  addRay(modelfluxAxis, 2.,155.5,2.5);

  CSpectrum spc = CSpectrum(std::move(spectralAxis),std::move(modelfluxAxis));

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
  BOOST_CHECK_CLOSE(gaussAmp, 0., precision);
  BOOST_CHECK_CLOSE(gaussPos, 0., precision);
  BOOST_CHECK_CLOSE(gaussWidth, 0., precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 0., precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 0., precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 0., precision);

  CGaussianFit::EStatus status = fitter.Compute( spc, TInt32Range( 0, 60 ) );
  BOOST_CHECK_EQUAL(status, CGaussianFit::EStatus::nStatus_Success);
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  fitter.GetResultsPolyCoeff0( coeff0 );
  BOOST_CHECK_CLOSE(gaussAmp, 1.5000000722300286, precision);
  BOOST_CHECK_CLOSE(gaussPos, 39.999999888990239, precision);
  BOOST_CHECK_CLOSE(gaussWidth, 4.0000003064126295, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 0.50559033514523466, precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 1.4705977499341338, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 1.7053687738517733, precision);
  BOOST_CHECK_CLOSE(coeff0, -9.4309808779202966e-08, 1e-06);

  status = fitter.Compute( spc, TInt32Range( 50, 50 ) );
  BOOST_CHECK_EQUAL(status,CGaussianFit::EStatus::nStatus_IllegalInput);

  status = fitter.Compute( spc, TInt32Range( 50, 52 ) );
  BOOST_CHECK_EQUAL(status,CGaussianFit::EStatus::nStatus_IllegalInput);

  status = fitter.Compute( spc, TInt32Range( 52, 60 ) );
  BOOST_CHECK_EQUAL(status,CGaussianFit::EStatus::nStatus_Success);

  status = fitter.Compute( spc, TInt32Range( 60, 100 ) );
  BOOST_CHECK_EQUAL(status, CGaussianFit::EStatus::nStatus_Success);
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  fitter.GetResultsPolyCoeff0( coeff0 );
  BOOST_CHECK_CLOSE(gaussAmp, 4.5000010751344064, precision);
  BOOST_CHECK_CLOSE(gaussPos, 80.000000178474721, precision);
  BOOST_CHECK_CLOSE(gaussWidth, 4.00000091535889, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 0.74797372318264976, precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 0.49396421482540065, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 0.79015532344560224, precision);
  BOOST_CHECK_CLOSE(coeff0, -1.2233871016444088e-06, 1e-06);

  status = fitter.Compute( spc, TInt32Range( 100, 140 ) );
  BOOST_CHECK_EQUAL(status, CGaussianFit::EStatus::nStatus_Success);
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  fitter.GetResultsPolyCoeff0( coeff0 );
  BOOST_CHECK_CLOSE(gaussAmp, -4.999997326237156, precision);
  BOOST_CHECK_CLOSE(gaussPos, 120.85000026571979, precision);
  BOOST_CHECK_CLOSE(gaussWidth, 3.1999982068897026, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 0.65213067912274125, precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 0.38948063181862164, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 0.52998315252751516, precision);
  BOOST_CHECK_CLOSE(coeff0, -3.1062808917571128e-06, 1e-06);

  status = fitter.Compute( spc, TInt32Range( 140, 170 ) );
  BOOST_CHECK_EQUAL(status, CGaussianFit::EStatus::nStatus_Success);
  fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
  fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );
  fitter.GetResultsPolyCoeff0( coeff0 );
  BOOST_CHECK_CLOSE(gaussAmp, 2.5000000545559202, precision);
  BOOST_CHECK_CLOSE(gaussPos, 155.50000007161785, precision);
  BOOST_CHECK_CLOSE(gaussWidth, 2.0000001160025973, precision);
  BOOST_CHECK_CLOSE(gaussAmpErr, 0.75192589874379412, precision);
  BOOST_CHECK_CLOSE(gaussPosErr, 0.60932628716452819, precision);
  BOOST_CHECK_CLOSE(gaussWidthErr, 0.76817236790330667, precision);
  BOOST_CHECK_CLOSE(coeff0, -8.6664245207542372e-08, 1e-06);
}

BOOST_AUTO_TEST_SUITE_END()

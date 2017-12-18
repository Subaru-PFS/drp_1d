#include <boost/test/unit_test.hpp>
#include <RedshiftLibrary/spectrum/spectralaxis.h>

using namespace NSEpic;

BOOST_AUTO_TEST_CASE(Constructor)
{
  // No element
  CSpectrumSpectralAxis n1Axis = CSpectrumSpectralAxis();
  BOOST_CHECK(n1Axis.GetSamplesCount()==0);
  // One element initialized to 0. linear scale
  CSpectrumSpectralAxis n21Axis = CSpectrumSpectralAxis(1, false);
  BOOST_CHECK(n21Axis.GetSamplesCount()==1);
  // One element initialized to 0. log scale
  CSpectrumSpectralAxis n22Axis = CSpectrumSpectralAxis(1, true);
  BOOST_CHECK(n22Axis.GetSamplesCount()==1);
  // One element initialized to one element array
  Float64 n3Array[] = {1.};
  CSpectrumSpectralAxis n31Axis = CSpectrumSpectralAxis(n3Array, 1, false);
  BOOST_CHECK(n31Axis.GetSamplesCount()==1);
  CSpectrumSpectralAxis n32Axis = CSpectrumSpectralAxis(n3Array, 1, true);
  BOOST_CHECK(n32Axis.GetSamplesCount()==1);
  // Clone from an other spectral axis redshift 0
  CSpectrumSpectralAxis n41Axis  = CSpectrumSpectralAxis(1, false);
  // CSpectrumSpectralAxis n42Axis  = CSpectrumSpectralAxis(n41Axis,0.,CSpectrumSpectralAxis::nShiftForward);
  // One element from 2 elements array
  Float64 n5Array[] = {1.,2.};
  CSpectrumSpectralAxis n5Axis = CSpectrumSpectralAxis(n5Array, 1, false);
  BOOST_CHECK(n5Axis.GetSamplesCount()==1);
  // Two elements from 1 element array
  Float64 n6Array[] = {1.};
  CSpectrumSpectralAxis n6Axis = CSpectrumSpectralAxis(n6Array, 2, false);
  BOOST_CHECK(n6Axis.GetSamplesCount()==2);
}

BOOST_AUTO_TEST_CASE(Operator)
{
  CSpectrumSpectralAxis n71Axis = CSpectrumSpectralAxis();
  Float64 n7Array[] = {1.};
  CSpectrumSpectralAxis n72Axis = CSpectrumSpectralAxis(n7Array, 1, false);
  n71Axis = n72Axis;
  BOOST_CHECK_CLOSE(n71Axis[0],n72Axis[0],1.e-12);
}

BOOST_AUTO_TEST_CASE(Resolution)
{
  Float64 n8Array[] = {1.,3.,4.,10.};
  CSpectrumSpectralAxis n8Axis = CSpectrumSpectralAxis(n8Array,4,false);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(),2.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetMeanResolution(),3.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(-1.),2.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(0.),2.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(1.),2.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(2.),2.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(3.),1.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(4.),6.,1.e-12);
  BOOST_CHECK_CLOSE(n8Axis.GetResolution(11.),6.,1.e-12);
}

BOOST_AUTO_TEST_CASE(ShiftByWaveLength)
{
  Float64 n9Array[] = {5.};
  CSpectrumSpectralAxis n91Axis = CSpectrumSpectralAxis(n9Array,1,false);
  n91Axis.ShiftByWaveLength(1.,CSpectrumSpectralAxis::nShiftForward);
  // BOOST_CHECK_CLOSE(n91Axis[0],6.,1.e-12);
  CSpectrumSpectralAxis n92Axis = CSpectrumSpectralAxis(n9Array,1,false);
  n92Axis.ShiftByWaveLength(1.,CSpectrumSpectralAxis::nShiftBackward);
  // BOOST_CHECK_CLOSE(n92Axis[0],4.,1.e-12);
  CSpectrumSpectralAxis n93Axis = CSpectrumSpectralAxis(n9Array,1,false);
  CSpectrumSpectralAxis n94Axis = CSpectrumSpectralAxis();
  CSpectrumSpectralAxis n95Axis = CSpectrumSpectralAxis();
  // n94Axis.ShiftByWaveLength(n93Axis,1.,CSpectrumSpectralAxis::nShiftForward);
  // n95Axis.ShiftByWaveLength(n93Axis,1.,CSpectrumSpectralAxis::nShiftBackward);
}

BOOST_AUTO_TEST_CASE(GetIndexAtWaveLength)
{
  Float64 n100Array[] = {0.,2.,3.,6.};
  CSpectrumSpectralAxis n100Axis = CSpectrumSpectralAxis(n100Array,4,false);
  BOOST_CHECK(n100Axis.GetIndexAtWaveLength(-1.)==0);
  BOOST_CHECK(n100Axis.GetIndexAtWaveLength(1.)==0);
  // BOOST_CHECK(n100Axis.GetIndexAtWaveLength(1.9)==1);
  // BOOST_CHECK(n100Axis.GetIndexAtWaveLength(2.1)==1);
  BOOST_CHECK(n100Axis.GetIndexAtWaveLength(7.)==3);
  BOOST_TEST_MESSAGE("index:"<<n100Axis.GetIndexAtWaveLength(-1.));
  BOOST_TEST_MESSAGE("index:"<<n100Axis.GetIndexAtWaveLength(1.));
  BOOST_TEST_MESSAGE("index:"<<n100Axis.GetIndexAtWaveLength(1.9));
  BOOST_TEST_MESSAGE("index:"<<n100Axis.GetIndexAtWaveLength(2.0));
  BOOST_TEST_MESSAGE("index:"<<n100Axis.GetIndexAtWaveLength(2.1));
  BOOST_TEST_MESSAGE("index:"<<n100Axis.GetIndexAtWaveLength(7.));
  TFloat64Range range1 = TFloat64Range(1.9,3.1);
  TInt32Range irange1  = n100Axis.GetIndexesAtWaveLengthRange(range1);
  BOOST_TEST_MESSAGE("index:"<<irange1.GetBegin()<<","<<irange1.GetEnd());
  BOOST_CHECK(irange1.GetBegin()==0);
  BOOST_CHECK(irange1.GetEnd()==3);
  TFloat64Range range2 = TFloat64Range(-2.0,-1.0);
  TInt32Range irange2  = n100Axis.GetIndexesAtWaveLengthRange(range2);
  BOOST_TEST_MESSAGE("index:"<<irange2.GetBegin()<<","<<irange2.GetEnd());
  BOOST_CHECK(irange2.GetBegin()==0);
  BOOST_CHECK(irange2.GetEnd()==0);
  TFloat64Range range3 = TFloat64Range(10.,20.);
  TInt32Range irange3  = n100Axis.GetIndexesAtWaveLengthRange(range3);
  BOOST_TEST_MESSAGE("index:"<<irange3.GetBegin()<<","<<irange3.GetEnd());
  BOOST_CHECK(irange3.GetBegin()==3);
  BOOST_CHECK(irange3.GetEnd()==3);
  TFloat64Range range4 = TFloat64Range(-1,2.2);
  TInt32Range irange4  = n100Axis.GetIndexesAtWaveLengthRange(range4);
  BOOST_TEST_MESSAGE("index:"<<irange4.GetBegin()<<","<<irange4.GetEnd());
  BOOST_CHECK(irange4.GetBegin()==0);
  BOOST_CHECK(irange4.GetEnd()==2);
  TFloat64Range range5 = TFloat64Range(2.2,10.);
  TInt32Range irange5  = n100Axis.GetIndexesAtWaveLengthRange(range5);
  BOOST_TEST_MESSAGE("index:"<<irange5.GetBegin()<<","<<irange5.GetEnd());
  BOOST_CHECK(irange5.GetBegin()==2);
  BOOST_CHECK(irange5.GetEnd()==3);
}

BOOST_AUTO_TEST_CASE(Scale)
{
  Float64 n111Array[] = {1.,1.};
  CSpectrumSpectralAxis n111Axis = CSpectrumSpectralAxis(n111Array,2,false);
  Float64 n111begin = n111Axis[0];
  BOOST_CHECK(n111Axis.IsInLinearScale()==true);
  BOOST_CHECK(n111Axis.IsInLogScale()==false);
  n111Axis.ConvertToLinearScale();
  BOOST_CHECK(n111Axis.IsInLinearScale()==true);
  BOOST_CHECK(n111Axis.IsInLogScale()==false);
  n111Axis.ConvertToLogScale();
  BOOST_CHECK(n111Axis.IsInLinearScale()==false);
  BOOST_CHECK(n111Axis.IsInLogScale()==true);
  BOOST_CHECK_CLOSE(n111Axis[0],0.,1e-12);

  CSpectrumSpectralAxis n112Axis = CSpectrumSpectralAxis(2,true);
  Float64 n112begin = n112Axis[0];
  BOOST_CHECK(n112Axis.IsInLinearScale()==false);
  BOOST_CHECK(n112Axis.IsInLogScale()==true);
  n112Axis.ConvertToLogScale();
  BOOST_CHECK(n112Axis.IsInLinearScale()==false);
  BOOST_CHECK(n112Axis.IsInLogScale()==true);
  n112Axis.ConvertToLinearScale();
  BOOST_CHECK(n112Axis.IsInLinearScale()==true);
  BOOST_CHECK(n112Axis.IsInLogScale()==false);
  BOOST_CHECK_CLOSE(n112Axis[0],1.,1e-12);
}

BOOST_AUTO_TEST_CASE(lambdaRange)
{
  Float64 n121Array[] = {0.,1.};
  CSpectrumSpectralAxis n121Axis = CSpectrumSpectralAxis(n121Array,2,false);
  TLambdaRange irange6 = n121Axis.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange6.GetBegin(),n121Array[0],1.e-12);
  BOOST_CHECK_CLOSE(irange6.GetEnd(),n121Array[1],1.e-12);
  CSpectrumSpectralAxis n122Axis = CSpectrumSpectralAxis(1,false);
  TLambdaRange irange7 = n122Axis.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange7.GetBegin(),0.,1.e-12);
  BOOST_CHECK_CLOSE(irange7.GetEnd(),0.,1.e-12);
  CSpectrumSpectralAxis n123Axis = CSpectrumSpectralAxis();
  TLambdaRange irange8 = n123Axis.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange8.GetBegin(),0.,1.e-12);
  BOOST_CHECK_CLOSE(irange8.GetEnd(),0.,1.e-12);
  Float64 n124Array[] = {0.};
  CSpectrumSpectralAxis n124Axis = CSpectrumSpectralAxis(n124Array,2,true);
  TLambdaRange irange9 = n124Axis.GetLambdaRange();
  BOOST_CHECK_CLOSE(irange9.GetBegin(),1.,1.e-12);
  BOOST_CHECK_CLOSE(irange9.GetEnd(),1.,1.e-12);
}

BOOST_AUTO_TEST_CASE(ClampLambdaRange)
{
  // Axis from [0. ... 1.]
  Float64 n131Array[] = {0.,1.};
  CSpectrumSpectralAxis n131Axis = CSpectrumSpectralAxis(n131Array,2,false);
  TFloat64Range crange;
  // Range [0.1 ... 0.9] inside axis range
  TFloat64Range irange10 = TFloat64Range(0.1,0.9);
  BOOST_CHECK(n131Axis.ClampLambdaRange(irange10, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(),0.1,1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(),0.9,1.e-12);
  // Range empty [0. ... 0.]
  TFloat64Range irange11 = TFloat64Range(0.,0.);
  BOOST_CHECK(~n131Axis.ClampLambdaRange(irange11, crange));
  // Axis empty
  CSpectrumSpectralAxis n132Axis = CSpectrumSpectralAxis(2,false);
  BOOST_CHECK(~n132Axis.ClampLambdaRange(irange10, crange));
  // Range [-1.0 ... 0.9] starts before axis [0. ... 1.]
  TFloat64Range irange12 = TFloat64Range(-1.0,0.9);
  BOOST_CHECK(n131Axis.ClampLambdaRange(irange12, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(),0,1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(),0.9,1.e-12);
  // Range [0.1 ... 1.1] ends after axis [0. ... 1.]
  TFloat64Range irange13 = TFloat64Range(0.1,1.1);
  BOOST_CHECK(n131Axis.ClampLambdaRange(irange13, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(),0.1,1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(),1.,1.e-12);
  // Range [-2. ... -1.] outside and before axis range
  TFloat64Range irange14 = TFloat64Range(-2.,-1.);
  BOOST_CHECK(n131Axis.ClampLambdaRange(irange14, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(),0.,1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(),-1.,1.e-12);
  // Range [-2. ... -1.] outside and after axis range
  TFloat64Range irange15 = TFloat64Range(2.,3.);
  BOOST_CHECK(n131Axis.ClampLambdaRange(irange15, crange));
  BOOST_CHECK_CLOSE(crange.GetBegin(),2.,1.e-12);
  BOOST_CHECK_CLOSE(crange.GetEnd(),1.,1.e-12);
}

  // CSpectrumSpectralAxis c42Axis  = CSpectrumSpectralAxis(n41Axis,0.,CSpectrumSpectralAxis::nShiftBackward);
  // BOOST_CHECK(c42Axis.GetSamplesCount()==1);
/*
  Float64 n1Array[] = {0.5};
  CSpectrumSpectralAxis n1Axis = CSpectrumSpectralAxis(n1Array, 1, false);
  BOOST_CHECK(n1Axis.GetSamples()[0] == n1Array[0]);

  Float64 n2Array[] = {0.,1.};
  CSpectrumSpectralAxis n2Axis = CSpectrumSpectralAxis(n2Array, 2, false);
  BOOST_CHECK(n2Axis.GetSamples()[0] == n2Array[0]);
  BOOST_CHECK(n2Axis.GetSamples()[1] == n2Array[1]);
*/

  //BOOST_TEST_MESSAGE("Resolution:"<<xAxis.GetResolution());
  //BOOST_TEST_MESSAGE("Mean resolution:"<<xAxis.GetMeanResolution());
  //BOOST_TEST_MESSAGE("Interpolate value:"<<yInt[0]);
  //BOOST_CHECK( 0 == 1 );

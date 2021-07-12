#include "RedshiftLibrary/spectrum/fluxaxis.h"

#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/mean.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
#include "RedshiftLibrary/spectrum/spectrum.h"

#include "RedshiftLibrary/log/log.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <numeric>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Spectrum)

bool correctMessage(const std::runtime_error& ex)
{
    BOOST_CHECK_EQUAL(ex.what(), std::string("Cannot rebin spectrum: target spectral axis is not included in the current spectral axis"));
    return true;
}
BOOST_AUTO_TEST_CASE(calcul)
{

  //--------------------//
  // constructor

  CSpectrumFluxAxis object_FluxAxis;
  BOOST_CHECK(object_FluxAxis.GetSamplesCount() == 0);

  UInt32 n = 10;
  CSpectrumFluxAxis *object_FluxAxis2 = new CSpectrumFluxAxis(n);
  BOOST_CHECK(object_FluxAxis2->GetSamplesCount() == n);

  Float64 Array1[10] = {1,2,3,4,5,6,7,8,9,10};
  CSpectrumFluxAxis *object_FluxAxis3 = new CSpectrumFluxAxis(Array1, 10);
  BOOST_CHECK(object_FluxAxis3->GetSamplesCount() == 10);
  BOOST_CHECK(accumulate(object_FluxAxis3->GetSamples(),
                         object_FluxAxis3->GetSamples() + 10, 0) == 55);
  BOOST_TEST_MESSAGE("index2:" << object_FluxAxis3->GetSamplesCount());

  Float64 Array2[10] = {2,4,6,8,10,12,14,16,18,20};
  CSpectrumFluxAxis *object_FluxAxis4 =
    new CSpectrumFluxAxis(Array1, 10, Array2, 10);
  Int32 s = object_FluxAxis4->GetSamplesCount();
  BOOST_CHECK(object_FluxAxis4->GetSamplesCount() == 10);

  //----------//
  // test rbin


  TFloat64Range *object_Range = new TFloat64Range(1., 5.);
  CSpectrumFluxAxis *sourceFluxAxis = new CSpectrumFluxAxis(Array1, n);
  CSpectrumSpectralAxis *sourceSpectralAxis =
    new CSpectrumSpectralAxis(n, false);
  CSpectrumSpectralAxis *targetSpectralAxis =
    new CSpectrumSpectralAxis(n, false);
  CSpectrumFluxAxis *rebinedFluxAxis = new CSpectrumFluxAxis(Array1, n);
  CSpectrumSpectralAxis *rebinedSpectralAxis =
    new CSpectrumSpectralAxis(n, false);
  CMask *rebinedMask = new CMask(n);

  // cas 1
  std::shared_ptr<CSpectrum> object_CSpectrum;
  object_CSpectrum = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis), (*sourceFluxAxis)));
  CSpectrum rebinnedSpectrum;
  bool resultcas1 = object_CSpectrum->Rebin((*object_Range),
						    (*targetSpectralAxis),
						    rebinnedSpectrum,
						    (*rebinedMask));
  BOOST_CHECK(resultcas1 == true);
  BOOST_TEST_MESSAGE("cas1:" << object_CSpectrum->Rebin((*object_Range),
								(*targetSpectralAxis),
								rebinnedSpectrum,
								(*rebinedMask)));


  // cas 3

  CSpectrumSpectralAxis *sourceSpectralAxis3 =
    new CSpectrumSpectralAxis(n, true);

  std::shared_ptr<CSpectrum> object_CSpectrum3;
  object_CSpectrum3 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis3), (*sourceFluxAxis)));
  CSpectrum rebinnedSpectrum3;
  bool resultcas3 = object_CSpectrum3->Rebin((*object_Range),
						    (*targetSpectralAxis),
						    rebinnedSpectrum3,
						    (*rebinedMask));
  BOOST_CHECK(resultcas3 == false);
  BOOST_TEST_MESSAGE("cas3:" << resultcas3);

  // cas 4

  TFloat64Range *object_Range4 = new TFloat64Range(3., 6.);
  TFloat64Range logIntersectedLambdaRange(log(object_Range4->GetBegin()),
                                          log(object_Range4->GetEnd()));
  TFloat64Range currentRange = logIntersectedLambdaRange;

  TFloat64List Array4 = {0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,
                      1.2,  1.3,  1.4,  1.45, 1.50, 1.55, 1.60,
                      1.65, 1.70, 1.75, 1.80, 1.85, 1.9};
  TFloat64List Array4b = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
  CSpectrumSpectralAxis *targetSpectralAxis4 =
    new CSpectrumSpectralAxis(Array4, true);
  CSpectrumSpectralAxis *sourceSpectralAxis4 =
    new CSpectrumSpectralAxis(Array4b, true);
  CMask *rebinedMask2 = new CMask();

  std::shared_ptr<CSpectrum> object_CSpectrum4;
  object_CSpectrum4 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis4), (*sourceFluxAxis)));
  CSpectrum rebinnedSpectrum4;
  BOOST_CHECK_EXCEPTION( object_CSpectrum4->Rebin((*object_Range4),(*targetSpectralAxis4),rebinnedSpectrum4,(*rebinedMask2)),
                        std::runtime_error, correctMessage);
  //BOOST_CHECK(resultcas4 == false);//we cannot interpolate outside spectrum range
  //BOOST_TEST_MESSAGE("cas4:" << resultcas4);

  //----------//
  // test rbin2

  // cas 1
  Float64 source = 1;
  const std::string opt_interp = "lin";

  std::shared_ptr<CSpectrum> object_CSpectrum5;
  object_CSpectrum5 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis), (*sourceFluxAxis)));
  CMask *rebinedMask3 = new CMask();
  CSpectrum rebinnedSpectrum5;
  bool resultRebincas1 = object_CSpectrum5->Rebin((*object_Range),
							   (*targetSpectralAxis),
							   rebinnedSpectrum5,
							   (*rebinedMask3),
							   opt_interp);
  BOOST_CHECK(resultRebincas1 == true);
  BOOST_TEST_MESSAGE("Rebin cas1:" << resultRebincas1);

  // cas 3
  CMask *rebinedMask5 = new CMask();
  std::shared_ptr<CSpectrum> object_CSpectrum7;
  object_CSpectrum7 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis3), (*sourceFluxAxis)));
  CSpectrum rebinnedSpectrum7;
  bool resultRebincas3 = object_CSpectrum7->Rebin(
							   (*object_Range), (*targetSpectralAxis), rebinnedSpectrum7, 
                 (*rebinedMask5), opt_interp);
  BOOST_CHECK(resultRebincas3 == false);
  BOOST_TEST_MESSAGE("Rebin cas3:" << resultRebincas3);

  // cas 4
    TFloat64List Array5b = {0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.1,
                      1.2,  1.3,  1.4,  1.45, 1.50, 1.55, 1.60,
                      1.65, 1.70, 1.75, 1.80, 1.85, 1.9};
  TFloat64List Array5 = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
  CSpectrumSpectralAxis *targetSpectralAxis5 =
    new CSpectrumSpectralAxis(Array5, true);
  CSpectrumSpectralAxis *sourceSpectralAxis5 =
    new CSpectrumSpectralAxis(Array5b, true);
  CSpectrumFluxAxis *sourceFluxAxis5 = new CSpectrumFluxAxis(Array1, n+n);

  std::shared_ptr<CSpectrum> object_CSpectrum8;
  object_CSpectrum8 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis5), (*sourceFluxAxis5)));
  CMask *rebinedMask6 = new CMask();
  CSpectrum rebinnedSpectrum8;
  bool resultRebincas4 = object_CSpectrum8->Rebin(
							   (*object_Range4), (*targetSpectralAxis5), 
                 rebinnedSpectrum8, (*rebinedMask6), opt_interp);
  BOOST_CHECK(resultRebincas4 == true);
  BOOST_TEST_MESSAGE("Rebin cas4:" << resultRebincas4);

  // cas 5
  const std::string opt_interp2 = "precomputedfinegrid";
  
  std::shared_ptr<CSpectrum> object_CSpectrum9;
  object_CSpectrum9 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis5), (*sourceFluxAxis5)));
  
  CMask *rebinedMask7 = new CMask();
  CSpectrum rebinnedSpectrum9;
  bool resultRebincas5 = object_CSpectrum9->Rebin(
							   (*object_Range4), (*targetSpectralAxis5), 
                 rebinnedSpectrum9, (*rebinedMask7), opt_interp2 );
  BOOST_CHECK(resultRebincas5 == true);
  BOOST_TEST_MESSAGE("Rebin cas5:" << resultRebincas5);

  CMask *rebinedMask8 = new CMask();
  std::shared_ptr<CSpectrum> object_CSpectrum10;
  object_CSpectrum10 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis5), (*sourceFluxAxis5)));
  std::shared_ptr<CSpectrum> object_CSpectrum_null = std::shared_ptr<CSpectrum>( new CSpectrum());
  CSpectrum rebinnedSpectrum10;
  bool resultRebincas5bis = object_CSpectrum10->Rebin(
							      (*object_Range4), (*targetSpectralAxis4), 
                    rebinnedSpectrum10, (*rebinedMask8), opt_interp2);
  BOOST_CHECK(resultRebincas5bis == true);
  BOOST_TEST_MESSAGE("Rebin cas5bis:" << resultRebincas5bis);
  // cas 6
  const std::string opt_interp3 = "spline";

  CMask *rebinedMask9 = new CMask();
  std::shared_ptr<CSpectrum> object_CSpectrum11;
  object_CSpectrum11 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis5), (*sourceFluxAxis5)));
  CSpectrum rebinnedSpectrum11;
  bool resultRebincas6 = object_CSpectrum11->Rebin(
							   (*object_Range4), (*targetSpectralAxis4), 
                 rebinnedSpectrum11, (*rebinedMask9), opt_interp3 );
  BOOST_CHECK(resultRebincas6 == true);
  BOOST_TEST_MESSAGE("Rebin cas6:" << resultRebincas6);

  // cas 7
  const std::string opt_interp4 = "ngp";
  CMask *rebinedMask10 = new CMask();
  std::shared_ptr<CSpectrum> object_CSpectrum12;
  object_CSpectrum12 = std::shared_ptr<CSpectrum>( new CSpectrum( (*sourceSpectralAxis5), (*sourceFluxAxis5)));
  CSpectrum rebinnedSpectrum12;
  bool resultRebincas7 = object_CSpectrum12->Rebin(
							   (*object_Range4), (*targetSpectralAxis4), 
                 rebinnedSpectrum12, (*rebinedMask10), opt_interp4);
  BOOST_CHECK(resultRebincas7 == true);
  BOOST_TEST_MESSAGE("Rebin cas7:" << resultRebincas7);

  //---------//
  // test ApplyMedianSmooth
  CSpectrumFluxAxis *object_CSpectrumFluxAxis =
    new CSpectrumFluxAxis(Array1, n);

  bool resultApplyMedianSmoothcas1 =
    object_CSpectrumFluxAxis->ApplyMedianSmooth(0);
  BOOST_CHECK(resultApplyMedianSmoothcas1 == false);
  BOOST_TEST_MESSAGE(
		     "resultApplyMedianSmooth cas1:" << resultApplyMedianSmoothcas1);

  BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis->GetSamplesCount"
                     << object_CSpectrumFluxAxis->GetSamplesCount());

  bool resultApplyMedianSmoothcas2 =
    object_CSpectrumFluxAxis->ApplyMedianSmooth(15);
  BOOST_CHECK(resultApplyMedianSmoothcas2 == false);
  BOOST_TEST_MESSAGE(
		     "resultApplyMedianSmooth cas2:" << resultApplyMedianSmoothcas2);

  bool resultApplyMedianSmoothcas3 =
    object_CSpectrumFluxAxis->ApplyMedianSmooth(5);
  BOOST_CHECK(resultApplyMedianSmoothcas3 == true);
  BOOST_TEST_MESSAGE(
		     "resultApplyMedianSmooth cas3:" << resultApplyMedianSmoothcas3);

  //--------------------//
  // ApplyMeanSmooth

  bool resultApplyMeanSmoothcas1 = object_CSpectrumFluxAxis->ApplyMeanSmooth(0);
  BOOST_CHECK(resultApplyMeanSmoothcas1 == false);
  BOOST_TEST_MESSAGE(
		     "resultApplyMeanSmooth cas1:" << resultApplyMeanSmoothcas1);

  BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis->GetSamplesCount"
                     << object_CSpectrumFluxAxis->GetSamplesCount());

  bool resultApplyMeanSmoothcas2 =
    object_CSpectrumFluxAxis->ApplyMeanSmooth(15);
  BOOST_CHECK(resultApplyMeanSmoothcas2 == false);
  BOOST_TEST_MESSAGE(
		     "resultApplyMeanSmooth cas2:" << resultApplyMeanSmoothcas2);

  bool resultApplyMeanSmoothcas3 = object_CSpectrumFluxAxis->ApplyMeanSmooth(5);
  BOOST_CHECK(resultApplyMeanSmoothcas3 == true);
  BOOST_TEST_MESSAGE(
		     "resultApplyMeanSmooth cas3:" << resultApplyMeanSmoothcas3);

  //--------------------//
  // test ComputeRMSDiff

  Float64 ArrayA[] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
  CSpectrumFluxAxis *object_FluxAxisA = new CSpectrumFluxAxis(ArrayA, 10);

  Float64 ArrayB[] = {2., 4., 6., 8., 10., 12., 14., 16., 18., 20.};
  CSpectrumFluxAxis *object_FluxAxisB = new CSpectrumFluxAxis(ArrayB, 10);

  Float64 resultComputeRMSDiff =
    object_FluxAxisA->ComputeRMSDiff((*object_FluxAxisB));
  BOOST_TEST_MESSAGE(
		     "object_FluxAxisA.GetSamples()[1]=" << object_FluxAxisA->GetSamples()[1]);
  BOOST_TEST_MESSAGE(
		     "object_FluxAxisB.GetSamples()[1]=" << object_FluxAxisB->GetSamples()[1]);

  Float64 er2 = 0.f;
  Float64 er = 0.f;

  Float64 weight = (Float64)n;
  for (int j = 0; j < 10; j++)
    {
      er2 += (object_FluxAxisA->GetSamples()[j] -
              object_FluxAxisB->GetSamples()[j]) *
	(object_FluxAxisA->GetSamples()[j] -
	 object_FluxAxisB->GetSamples()[j]) /
	weight;
    }
  er = sqrt(er2);

  BOOST_CHECK_CLOSE(resultComputeRMSDiff, er, 1.e-12);

  //--------------------//
  // test Subtract

  Bool resultSubtract = object_FluxAxisA->Subtract((*object_FluxAxisB));

  int indice = 0;
  bool sub;

  for (UInt32 i = 0; i < 10; i++)
    {
      object_FluxAxisA->GetSamples()[i] =
	object_FluxAxisA->GetSamples()[i] - object_FluxAxisB->GetSamples()[i];

      indice++;
    }

  if (indice == 10)
    {

      sub = true;
    }
  else
    {

      sub = false;
    }

  BOOST_CHECK(resultSubtract == sub);
  BOOST_TEST_MESSAGE("resultSubtract=" << resultSubtract << ", sub=" << sub);

  //--------------------//
  // test Invert

  Bool resultInvert = object_FluxAxisA->Invert();

  int indice2 = 0;
  bool inv;

  for (UInt32 i = 0; i < 10; i++)
    {
      object_FluxAxisA->GetSamples()[i] = -object_FluxAxisA->GetSamples()[i];

      indice2++;
    }

  if (indice2 == 10)
    {

      inv = true;
    }
  else
    {

      inv = false;
    }

  BOOST_CHECK(resultInvert == inv);
  BOOST_TEST_MESSAGE("resultInvert=" << resultInvert << ", inv=" << inv);

  //--------------------//
  // test ComputeMeanAndSDevWithoutError

  CMask *Mask = new CMask(10);

  Float64 mean = 1.;
  Float64 sdev = 1.;
  const CSpectrumNoiseAxis error(1);
  const CSpectrumNoiseAxis empty_error;

  object_FluxAxisA->GetError() = error;
  Bool resultComputeMeanAndSDev_cas1 =
    object_FluxAxisA->ComputeMeanAndSDev((*Mask), mean, sdev);
  BOOST_CHECK(resultComputeMeanAndSDev_cas1 == false);

  for (int i=0;i<10;i++) {
    (*Mask)[i] = 0;
  }

  object_FluxAxisA->GetError() = empty_error;
  Bool resultComputeMeanAndSDev_cas2 =
    object_FluxAxisA->ComputeMeanAndSDev((*Mask), mean, sdev);
  BOOST_CHECK(resultComputeMeanAndSDev_cas2 == false);
}

BOOST_AUTO_TEST_CASE(RebinVarianceWeighted)
{
  // test RebinVarianceWeighted
  TFloat64List lambdas = {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000};
  CSpectrumFluxAxis sourceFluxAxis(10);
  TFloat64List Array = {1000, 2500, 2900, 3890, 4690, 5500, 6800, 7001, 8033, 10000};
  CSpectrumSpectralAxis sourceSpectralAxis(Array, false);
  CSpectrumSpectralAxis bogus_sourceSpectralAxis(4, false);
 
  CSpectrumSpectralAxis targetSpectralAxis(lambdas);
  TFloat64Range currentRange = targetSpectralAxis.GetLambdaRange();
  
  const std::string opt_interp = "lin";
  std::string errorRebinMethod = "rebinVariance";

  std::shared_ptr<CSpectrum> object_CSpectrum13;
  object_CSpectrum13 = std::shared_ptr<CSpectrum>( new CSpectrum( (sourceSpectralAxis), (sourceFluxAxis)));
  CMask *rebinedMask13 = new CMask();
  CSpectrum rebinnedSpectrum13;
  bool resultRebinvaria1 = object_CSpectrum13->Rebin(currentRange,
							   (targetSpectralAxis),
							   rebinnedSpectrum13,
							   (*rebinedMask13),
							   opt_interp, errorRebinMethod);
  BOOST_CHECK(resultRebinvaria1 == true);

}


BOOST_AUTO_TEST_SUITE_END()

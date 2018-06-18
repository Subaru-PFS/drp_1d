#include <RedshiftLibrary/spectrum/fluxaxis.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/common/median.h>
#include <RedshiftLibrary/common/mean.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>

#include <math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <RedshiftLibrary/log/log.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Spectrum)

BOOST_AUTO_TEST_CASE(calcul)
{

//--------------------//
//constructor

CSpectrumFluxAxis object_FluxAxis;
BOOST_CHECK(object_FluxAxis.GetSamplesCount()==0);

UInt32 n=10;
CSpectrumFluxAxis* object_FluxAxis2=new CSpectrumFluxAxis(n);
BOOST_CHECK(object_FluxAxis2->GetSamplesCount()==n);

Float64 Array1[] = {1.};
CSpectrumFluxAxis* object_FluxAxis3=new CSpectrumFluxAxis(Array1,n);
BOOST_CHECK(object_FluxAxis3->GetSamplesCount()==n);
BOOST_TEST_MESSAGE("index2:"<<object_FluxAxis3->GetSamplesCount());

//----------//
//test rbin
CSpectrumFluxAxis* object_CSpectrumFluxAxis=new CSpectrumFluxAxis(Array1,n);

TFloat64Range* object_Range=new TFloat64Range(1.,5.);
CSpectrumFluxAxis* sourceFluxAxis=new CSpectrumFluxAxis(Array1,n);
CSpectrumSpectralAxis* sourceSpectralAxis = new CSpectrumSpectralAxis(n, false);
CSpectrumSpectralAxis* targetSpectralAxis = new CSpectrumSpectralAxis(n, false);
CSpectrumFluxAxis* rebinedFluxAxis=new CSpectrumFluxAxis(Array1,n);
CSpectrumSpectralAxis* rebinedSpectralAxis = new CSpectrumSpectralAxis(n, false);
CMask* rebinedMask=new CMask;

//cas 1

bool resultcas1=object_CSpectrumFluxAxis->Rebin((*object_Range),(*sourceFluxAxis),(*sourceSpectralAxis),(*targetSpectralAxis),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask));
BOOST_CHECK(resultcas1==true);
BOOST_TEST_MESSAGE("cas1:"<<object_CSpectrumFluxAxis->Rebin((*object_Range),(*sourceFluxAxis),(*sourceSpectralAxis),(*targetSpectralAxis),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask)));

//cas 2

CSpectrumSpectralAxis* sourceSpectralAxis2 = new CSpectrumSpectralAxis(n+1, false);

bool resultcas2=object_CSpectrumFluxAxis->Rebin((*object_Range),(*sourceFluxAxis),(*sourceSpectralAxis2),(*targetSpectralAxis),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask));
BOOST_CHECK(resultcas2==false);
BOOST_TEST_MESSAGE("cas2:"<<resultcas2);

//cas 3

CSpectrumSpectralAxis* sourceSpectralAxis3 = new CSpectrumSpectralAxis(n, true);

bool resultcas3=object_CSpectrumFluxAxis->Rebin((*object_Range),(*sourceFluxAxis),(*sourceSpectralAxis3),(*targetSpectralAxis),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask));
BOOST_CHECK(resultcas3==false);
BOOST_TEST_MESSAGE("cas3:"<<resultcas3);

//cas 4

TFloat64Range* object_Range4=new TFloat64Range(3.,6.);
TFloat64Range logIntersectedLambdaRange( log( object_Range4->GetBegin() ), log( object_Range4->GetEnd() ) );
TFloat64Range currentRange = logIntersectedLambdaRange;

Float64 Array4[] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.45,1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.9};
Float64 Array4b[] = {1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9};
CSpectrumSpectralAxis* targetSpectralAxis4 = new CSpectrumSpectralAxis(Array4,n+n, true);
CSpectrumSpectralAxis* sourceSpectralAxis4 = new CSpectrumSpectralAxis(Array4b,n, true);

bool resultcas4=object_CSpectrumFluxAxis->Rebin((*object_Range4),(*sourceFluxAxis),(*sourceSpectralAxis4),(*targetSpectralAxis4),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask));
BOOST_CHECK(resultcas4==true);
BOOST_TEST_MESSAGE("cas4:"<<resultcas4);


//----------//
//test rbin2

//cas 1

const Float64* Buffer=sourceSpectralAxis->GetSamples();
Float64 source=1;
const std::string opt_interp="lin";

bool resultrebin2cas1=object_CSpectrumFluxAxis->Rebin2((*object_Range),(*sourceFluxAxis),Buffer,source,(*sourceSpectralAxis),(*targetSpectralAxis),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask),opt_interp);
BOOST_CHECK(resultrebin2cas1==true);
BOOST_TEST_MESSAGE("rebin2 cas1:"<<resultrebin2cas1);

//cas 2

bool resultrebin2cas2=object_CSpectrumFluxAxis->Rebin2((*object_Range),(*sourceFluxAxis),Buffer,source,(*sourceSpectralAxis2),(*targetSpectralAxis),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask),opt_interp);
BOOST_CHECK(resultrebin2cas2==false);
BOOST_TEST_MESSAGE("rebin2 cas2:"<<resultrebin2cas2);

//cas 3

bool resultrebin2cas3=object_CSpectrumFluxAxis->Rebin2((*object_Range),(*sourceFluxAxis),Buffer,source,(*sourceSpectralAxis3),(*targetSpectralAxis),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask),opt_interp);
BOOST_CHECK(resultrebin2cas3==false);
BOOST_TEST_MESSAGE("rebin2 cas3:"<<resultrebin2cas3);

//cas 4

bool resultrebin2cas4=object_CSpectrumFluxAxis->Rebin2((*object_Range4),(*sourceFluxAxis),Buffer,source,(*sourceSpectralAxis4),(*targetSpectralAxis4),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask),opt_interp);
BOOST_CHECK(resultrebin2cas4==true);
BOOST_TEST_MESSAGE("rebin2 cas4:"<<resultrebin2cas4);

//cas 5
const std::string opt_interp2="precomputedfinegrid";

bool resultrebin2cas5=object_CSpectrumFluxAxis->Rebin2((*object_Range4),(*sourceFluxAxis),Buffer,source,(*sourceSpectralAxis4),(*targetSpectralAxis4),(*rebinedFluxAxis),(*rebinedSpectralAxis),(*rebinedMask),opt_interp2);
BOOST_CHECK(resultrebin2cas5==true);
BOOST_TEST_MESSAGE("rebin2 cas5:"<<resultrebin2cas5);


//---------//
//test ApplyMedianSmooth

bool resultApplyMedianSmoothcas1=object_CSpectrumFluxAxis->ApplyMedianSmooth(0);
BOOST_CHECK(resultApplyMedianSmoothcas1==false);
BOOST_TEST_MESSAGE("resultApplyMedianSmooth cas1:"<<resultApplyMedianSmoothcas1);

BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis->GetSamplesCount"<<object_CSpectrumFluxAxis->GetSamplesCount());

bool resultApplyMedianSmoothcas2=object_CSpectrumFluxAxis->ApplyMedianSmooth(15);
BOOST_CHECK(resultApplyMedianSmoothcas2==false);
BOOST_TEST_MESSAGE("resultApplyMedianSmooth cas2:"<<resultApplyMedianSmoothcas2);

bool resultApplyMedianSmoothcas3=object_CSpectrumFluxAxis->ApplyMedianSmooth(5);
BOOST_CHECK(resultApplyMedianSmoothcas3==true);
BOOST_TEST_MESSAGE("resultApplyMedianSmooth cas3:"<<resultApplyMedianSmoothcas3);

//--------------------//
//ApplyMeanSmooth

bool resultApplyMeanSmoothcas1=object_CSpectrumFluxAxis->ApplyMeanSmooth(0);
BOOST_CHECK(resultApplyMeanSmoothcas1==false);
BOOST_TEST_MESSAGE("resultApplyMeanSmooth cas1:"<<resultApplyMeanSmoothcas1);

BOOST_TEST_MESSAGE("object_CSpectrumFluxAxis->GetSamplesCount"<<object_CSpectrumFluxAxis->GetSamplesCount());

bool resultApplyMeanSmoothcas2=object_CSpectrumFluxAxis->ApplyMeanSmooth(15);
BOOST_CHECK(resultApplyMeanSmoothcas2==false);
BOOST_TEST_MESSAGE("resultApplyMeanSmooth cas2:"<<resultApplyMeanSmoothcas2);

bool resultApplyMeanSmoothcas3=object_CSpectrumFluxAxis->ApplyMeanSmooth(5);
BOOST_CHECK(resultApplyMeanSmoothcas3==true);
BOOST_TEST_MESSAGE("resultApplyMeanSmooth cas3:"<<resultApplyMeanSmoothcas3);

//--------------------//
// test ComputeRMSDiff

Float64 ArrayA[] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
CSpectrumFluxAxis* object_FluxAxisA=new CSpectrumFluxAxis(ArrayA,10);

Float64 ArrayB[] = {2.,4.,6.,8.,10.,12.,14.,16.,18.,20.};
CSpectrumFluxAxis* object_FluxAxisB=new CSpectrumFluxAxis(ArrayB,10);

Float64 resultComputeRMSDiff=object_FluxAxisA->ComputeRMSDiff((*object_FluxAxisB));
BOOST_TEST_MESSAGE("object_FluxAxisA.GetSamples()[1]="<<object_FluxAxisA->GetSamples()[1]);
BOOST_TEST_MESSAGE("object_FluxAxisB.GetSamples()[1]="<<object_FluxAxisB->GetSamples()[1]);

Float64 er2 = 0.f;
Float64 er = 0.f;

Float64 weight = (Float64)n;
for (int j=0;j<10;j++)
{
    er2 += (object_FluxAxisA->GetSamples()[j]-object_FluxAxisB->GetSamples()[j])*(object_FluxAxisA->GetSamples()[j]-object_FluxAxisB->GetSamples()[j]) / weight;
}
er = sqrt(er2);

BOOST_CHECK_CLOSE(resultComputeRMSDiff,er,1.e-12);

//--------------------//
// test Subtract

Bool resultSubtract=object_FluxAxisA->Subtract((*object_FluxAxisB));

int indice=0;
bool sub;

for( UInt32 i=0; i<10; i++ )
{
    object_FluxAxisA->GetSamples()[i] = object_FluxAxisA->GetSamples()[i]-object_FluxAxisB->GetSamples()[i];

    indice++;
}

if(indice==10)
{

  sub=true;

}
else
{

  sub=false;

}

BOOST_CHECK(resultSubtract==sub);
BOOST_TEST_MESSAGE("resultSubtract="<<resultSubtract<<", sub="<<sub);

//--------------------//
// test Invert

Bool resultInvert=object_FluxAxisA->Invert();

int indice2=0;
bool inv;

for( UInt32 i=0; i<10; i++ )
{
    object_FluxAxisA->GetSamples()[i] = -object_FluxAxisA->GetSamples()[i];

    indice2++;
}

if(indice2==10)
{

  inv=true;

}
else
{

  inv=false;

}

BOOST_CHECK(resultInvert==inv);
BOOST_TEST_MESSAGE("resultInvert="<<resultInvert<<", inv="<<inv);

//--------------------//
//test ComputeMeanAndSDevWithoutError

CMask* Mask=new CMask(10);

// for (int i=0;i<10;i++)
// {
// Mask[i] = 0;
// }

Float64 mean=1.;
Float64 sdev=1.;
const TFloat64List error { 1. };

Bool resultComputeMeanAndSDev_cas1=object_FluxAxisA->ComputeMeanAndSDev((*Mask),mean,sdev,error);
BOOST_CHECK(resultComputeMeanAndSDev_cas1==false);


//--------------------//

}

BOOST_AUTO_TEST_SUITE_END()

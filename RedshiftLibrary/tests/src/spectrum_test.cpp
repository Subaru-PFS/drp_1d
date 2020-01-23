#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>


#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/debug/assert.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <limits>
#include <vector>

#include <gsl/gsl_fit.h>

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/common/datatypes.h>

#include <RedshiftLibrary/spectrum/io/fitsreader.h>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;


BOOST_AUTO_TEST_SUITE(Spectrum)

BOOST_AUTO_TEST_CASE(name)
{

    CSpectrum object_CSpectrum;

    const char* test_name="toto";
    object_CSpectrum.SetName(test_name);
    BOOST_CHECK(object_CSpectrum.GetName()==test_name);
    BOOST_TEST_MESSAGE("GetName OK");

}

BOOST_AUTO_TEST_CASE(path)
{

    CSpectrum object_CSpectrum;

    const char* test_path="chemin";
    object_CSpectrum.SetFullPath(test_path);
    BOOST_CHECK(object_CSpectrum.GetFullPath()==test_path);
    BOOST_TEST_MESSAGE("GetFullPath OK");

}

BOOST_AUTO_TEST_CASE(decompScales)
{

    CSpectrum object_CSpectrum;

    Int32 decompScales=10;
    object_CSpectrum.SetDecompScales(decompScales);
    BOOST_CHECK(object_CSpectrum.GetDecompScales()==decompScales);
    BOOST_TEST_MESSAGE("GetDecompScales OK");

}
BOOST_AUTO_TEST_CASE(invert)
{

    CSpectrum object_CSpectrum;

    BOOST_CHECK(object_CSpectrum.InvertFlux() == true);
    BOOST_TEST_MESSAGE("InvertFlux OK");

}

BOOST_AUTO_TEST_CASE(Calcul)
{

    //--------------------//
    //constructor

    CSpectrum object_CSpectrum;
    BOOST_TEST_MESSAGE("index:"<<object_CSpectrum.GetSampleCount());
    BOOST_CHECK(object_CSpectrum.GetSampleCount()==0);

    CSpectrumFluxAxis m_FluxAxis;
    CSpectrumSpectralAxis m_SpectralAxis;

    int nbmin=0;
    int nbmax=11;

    CSpectrumSpectralAxis *_SpectralAxis = new CSpectrumSpectralAxis(nbmax, false);
    CSpectrumFluxAxis *_FluxAxis = new CSpectrumFluxAxis(nbmax);

    for (int i=nbmin; i<nbmax;i++)
    {
        (*_SpectralAxis)[i] = i+1;

        if (i<(nbmax-nbmin)/2)
        {
            (*_FluxAxis)[i] = 0.0;
            (*_FluxAxis).GetError()[i] = 0.0;
        }
        else if (i==7)
        {
            (*_FluxAxis)[i] = std::nan("1");
            (*_FluxAxis).GetError()[i] = std::nan("2");
        }
        else if (i==9)
        {
            (*_FluxAxis)[i] = std::numeric_limits<double>::infinity();
            (*_FluxAxis).GetError()[i] = std::numeric_limits<double>::infinity();
        }
        else
        {
            (*_FluxAxis)[i] = i+2;
            (*_FluxAxis).GetError()[i] = 1e-12;
        }

        BOOST_TEST_MESSAGE("(*_SpectralAxis)[i]:"<<(*_SpectralAxis)[i]);
    }

    m_SpectralAxis = *_SpectralAxis;
    m_FluxAxis = *_FluxAxis;

    BOOST_TEST_MESSAGE("index1:"<<m_SpectralAxis.GetSamplesCount());
    BOOST_TEST_MESSAGE("index2:"<<m_FluxAxis.GetSamplesCount());

    object_CSpectrum.GetSpectralAxis()=m_SpectralAxis;
    object_CSpectrum.GetFluxAxis()=m_FluxAxis;

    BOOST_TEST_MESSAGE("index21:"<<object_CSpectrum.GetSampleCount());
    BOOST_CHECK(object_CSpectrum.GetSampleCount()==nbmax);

    CSpectrum object_CSpectrum2;

    object_CSpectrum2=object_CSpectrum;
    BOOST_CHECK(object_CSpectrum2.GetSampleCount()==nbmax);

    BOOST_TEST_MESSAGE("index22:"<<object_CSpectrum2.GetSampleCount());
    BOOST_TEST_MESSAGE("index23:"<<object_CSpectrum2.GetFluxAxis()[0]);
    BOOST_TEST_MESSAGE("index23:"<<object_CSpectrum2.GetSpectralAxis()[0]);

    TFloat64List mask;
    mask.push_back(1);
    BOOST_TEST_MESSAGE("index33:"<<mask.size());
    BOOST_TEST_MESSAGE("index33:"<<mask[0]);
    BOOST_TEST_MESSAGE("index34:"<<mask[0]);

    CSpectrum* object_CSpectrum3=new CSpectrum(object_CSpectrum2,mask);

    BOOST_CHECK(object_CSpectrum3->GetSampleCount()==1);
    BOOST_CHECK_CLOSE(object_CSpectrum3->GetFluxAxis()[0],object_CSpectrum2.GetFluxAxis()[0],1e-12);
    BOOST_CHECK_CLOSE(object_CSpectrum3->GetSpectralAxis()[0],object_CSpectrum2.GetSpectralAxis()[0],1e-12);


    BOOST_TEST_MESSAGE("index43:"<<object_CSpectrum3->GetFluxAxis()[0]);

    BOOST_TEST_MESSAGE("test constructeur OK");

    //--------------------//
    //test GetMeanAndStdFluxInRange

    TFloat64Range* object_CRange=new TFloat64Range(1.,5.);
    TFloat64Range* object_CRangeb=new TFloat64Range(-1.,5.);
    TFloat64Range* object_CRangec=new TFloat64Range(1.,15.);

    BOOST_TEST_MESSAGE("index51:"<<object_CRange->GetBegin());
    BOOST_TEST_MESSAGE("index52:"<<object_CRange->GetEnd());

    Float64 mean=10.0;
    Float64 std=10.0;

    bool result1=object_CSpectrum2.GetMeanAndStdFluxInRange((*object_CRange),mean,std);
    bool result1b=object_CSpectrum2.GetMeanAndStdFluxInRange((*object_CRangeb),mean,std);
    bool result1c=object_CSpectrum2.GetMeanAndStdFluxInRange((*object_CRangec),mean,std);

    BOOST_TEST_MESSAGE("index53:"<<result1);
    BOOST_TEST_MESSAGE("index53b:"<<result1b);
    BOOST_TEST_MESSAGE("index53c:"<<result1c);
    BOOST_TEST_MESSAGE("index54:"<<object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetBegin());
    BOOST_TEST_MESSAGE("index55:"<<object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetEnd());

    BOOST_CHECK( result1 == true );
    BOOST_CHECK( result1b == false );
    BOOST_CHECK( result1c == false );

    //--------------------//
    //test GetLinearRegInRange

    bool result11=object_CSpectrum2.GetLinearRegInRange((*object_CRange),mean,std);
    bool result11b=object_CSpectrum2.GetLinearRegInRange((*object_CRangeb),mean,std);
    bool result11c=object_CSpectrum2.GetLinearRegInRange((*object_CRangec),mean,std);

    BOOST_CHECK( result11 == true );
    BOOST_CHECK( result11b == false );
    BOOST_CHECK( result11c == false );

    //--------------------//
    //test IsFluxValid

    BOOST_CHECK(object_CSpectrum2.IsFluxValid(1,11)==false);//cas dans tout l'intervalle
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(1,5)==false);//cas dans l'intervalle 1 à 5 avec 0.0
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(1,6)==true);//cas dans l'intervalle 1 à 6
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(6,11)==false);//cas dans l'intervalle 6 à 11 avec nan et inf
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(6,7)==true);//cas dans l'intervalle 6 à 7
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(7,9)==false);//cas dans l'intervalle 7 à 9 avec nan
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(9,9)==true);//cas où l'intervalle est un point
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(9,11)==false);//cas dans l'intervalle 9 à 11 avec inf
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(10,10)==false);//cas où l'intervalle est un point inf
    BOOST_CHECK(object_CSpectrum2.IsFluxValid(11,14)==false);//cas où l'intervalle est à l'extérieur

    //--------------------//
    //test IsNoiseValid

    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(1,11)==false);//cas dans tout l'intervalle
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(1,5)==false);//cas dans l'intervalle 1 à 5 avec 0.0
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(1,6)==false);//cas dans l'intervalle 1 à 6
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(6,11)==false);//cas dans l'intervalle 6 à 11 avec nan et inf
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(6,7)==true);//cas dans l'intervalle 6 à 7
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(7,9)==false);//cas dans l'intervalle 7 à 9 avec nan
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(9,9)==true);//cas où l'intervalle est un point
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(9,11)==false);//cas dans l'intervalle 9 à 11 avec inf
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(10,10)==false);//cas où l'intervalle est un point inf
    BOOST_CHECK(object_CSpectrum2.IsNoiseValid(11,14)==false);//cas où l'intervalle est à l'extérieur

    //--------------------//
    //test correctSpectrum

    BOOST_CHECK(object_CSpectrum2.correctSpectrum(1,5)==false);//cas dans l'intervalle 1 à 5 avec 0.0
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(1,6)==true);//cas dans l'intervalle 1 à 6
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(6,7)==false);//cas dans l'intervalle 6 à 7
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(7,9)==true);//cas dans l'intervalle 7 à 9 avec nan
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(9,9)==false);//cas où l'intervalle est un point
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(9,11)==true);//cas dans l'intervalle 9 à 11 avec inf
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(10,10)==false);//cas où l'intervalle est un point inf
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(6,11)==false);//cas dans l'intervalle 6 à 11 avec nan et inf
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(1,11)==false);//cas dans tout l'intervalle
    BOOST_CHECK(object_CSpectrum2.correctSpectrum(11,14)==false);//cas où l'intervalle est à l'extérieur

    //--------------------//
    //test GetLambdaRange

    Float64 intervalle=object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetLength();

    BOOST_CHECK_CLOSE(intervalle,object_CSpectrum2.GetLambdaRange().GetLength(),1e-12);

    BOOST_TEST_MESSAGE("result9A:"<<object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetLength());
    BOOST_TEST_MESSAGE("index55A:"<<object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetBegin());
    BOOST_TEST_MESSAGE("index55B:"<<object_CSpectrum2.GetSpectralAxis().GetLambdaRange().GetEnd());

    //--------------------//
    //test ConvertToLinearScale

    bool result2=object_CSpectrum2.ConvertToLinearScale();
    bool result3=object_CSpectrum2.GetSpectralAxis().ConvertToLinearScale();

    BOOST_CHECK( result2 == result3 );

    BOOST_TEST_MESSAGE("result2:"<<result2);
    BOOST_TEST_MESSAGE("result3:"<<result3);


    //--------------------//
    ///test ConvertToLogScale

    bool result4=object_CSpectrum2.ConvertToLogScale();
    bool result5=object_CSpectrum2.GetSpectralAxis().ConvertToLogScale();

    BOOST_CHECK( result4 == result5 );

    BOOST_TEST_MESSAGE("result4:"<<result4);
    BOOST_TEST_MESSAGE("result5:"<<result5);


    //--------------------//
    //test GetMeanResolution

    Float64 result6=object_CSpectrum2.GetMeanResolution();
    Float64 result7=object_CSpectrum2.GetSpectralAxis().GetMeanResolution();

    BOOST_CHECK_CLOSE(result6,result7,1e-12);

    BOOST_TEST_MESSAGE("result6:"<<result6);
    BOOST_TEST_MESSAGE("result7:"<<result7);


    //--------------------//
    //test GetResolution

    Float64 result8=object_CSpectrum2.GetResolution();
    Float64 result9=object_CSpectrum2.GetSpectralAxis().GetResolution();

    BOOST_CHECK_CLOSE(result8,result9,1e-12);

    BOOST_TEST_MESSAGE("result8:"<<result8);
    BOOST_TEST_MESSAGE("result9:"<<result9);



    //--------------------//
    // //test removeContinuum

    CContinuumIrregularSamplingMedian remover2;

    BOOST_CHECK(object_CSpectrum2.RemoveContinuum(remover2)==true);
    BOOST_TEST_MESSAGE("test Remove:"<<object_CSpectrum2.RemoveContinuum(remover2));

    delete _SpectralAxis;
    delete _FluxAxis;
    delete object_CSpectrum3;
    delete object_CRange;
    delete object_CRangeb;
    delete object_CRangec;

}

BOOST_AUTO_TEST_SUITE_END()

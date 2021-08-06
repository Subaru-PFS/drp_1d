#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/LSFFactory.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
#include "RedshiftLibrary/spectrum/LSFVariableWidth.h"

#include "RedshiftLibrary/debug/assert.h"

#include <numeric>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <limits>
#include <vector>

#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/datatypes.h"
#include <boost/test/unit_test.hpp>

using namespace NSEpic;


BOOST_AUTO_TEST_SUITE(LSF)

bool correctMessage(const std::exception& ex)
{
    BOOST_CHECK_EQUAL(ex.what(), std::string("Size do not match "));
    return true;
}

BOOST_AUTO_TEST_CASE(LSF_ConstantWidth)
{

    CSpectrumSpectralAxis SpectralAxis;
    CSpectrumFluxAxis FluxAxis;
    std::string lsfType = "GaussianConstantWidth";
    Float64 width = 1.09;
    TScopeStack scopeStack;
    std::shared_ptr<CParameterStore> store = std::make_shared<CParameterStore>(scopeStack);
    store->Set( "LSF.width", width );
    std::shared_ptr<TLSFArguments> args = std::make_shared<TLSFGaussianConstantWidthArgs>(store);
    std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

    //Test constructor with spectralAxis, fluxAxis and LSF
    CSpectrum object_CSpectrum = CSpectrum(SpectralAxis, FluxAxis, LSF);
    Float64 lambda = 7000.;
    BOOST_CHECK(object_CSpectrum.GetLSF()->IsValid() == true);
    BOOST_CHECK(object_CSpectrum.GetLSF()->GetWidth(lambda) == 1.09);

    //Test assignment copy constructor
    CSpectrum object_CSpectrum1 = object_CSpectrum;

    BOOST_CHECK(object_CSpectrum1.GetLSF()->IsValid() == true);
    BOOST_CHECK(object_CSpectrum1.GetLSF()->GetWidth(lambda) == 1.09);

    //Test copy constructor
    CSpectrum object_CSpectrum1_bis(object_CSpectrum1);

    BOOST_CHECK(object_CSpectrum1_bis.GetLSF()->IsValid() == true);
    BOOST_CHECK(object_CSpectrum1_bis.GetLSF()->GetWidth(lambda) == 1.09);
    
    
    //Test constructor with spectralAxis and fluxAxis
    CSpectrum object_CSpectrum2(SpectralAxis, FluxAxis);
    BOOST_CHECK(object_CSpectrum2.GetLSF() == nullptr);
    /*
    //Test default constructor
    CSpectrum object_CSpectrum3;
    BOOST_CHECK(object_CSpectrum3.GetLSF() == nullptr);
    
    object_CSpectrum3.SetLSF(LSF);
    BOOST_CHECK(object_CSpectrum3.GetLSF()->IsValid() == true);
    BOOST_CHECK(object_CSpectrum3.GetLSF()->GetWidth() == 1.09);
    /*object_CSpectrum3.GetLSF()->SetWidth(2.04e-60);
    BOOST_CHECK(object_CSpectrum3.GetLSF()->GetWidth() == 2.04e-60); 
    object_CSpectrum3.GetLSF()->SetWidth(0.0);
    BOOST_CHECK(object_CSpectrum3.GetLSF()->IsValid() == false);
    BOOST_CHECK(object_CSpectrum3.GetLSF()->GetWidth() == 0.0);
    object_CSpectrum3.GetLSF()->SetWidth(DBL_MAX);
    BOOST_CHECK(object_CSpectrum3.GetLSF()->IsValid() == true);
    BOOST_CHECK(object_CSpectrum3.GetLSF()->GetWidth() == DBL_MAX);
    BOOST_TEST_MESSAGE("LSF OK");
    */
}

BOOST_AUTO_TEST_CASE(LSF_VariableWidth)
{/*
    TScopeStack scopeStack;
    std::shared_ptr<CParameterStore> store = std::make_shared<CParameterStore>(scopeStack);
*/
    Int32 n = 10;
    TFloat64List spcAxis(n);
    for(Int32 i =0; i<n; i++)
        spcAxis[i]= 100+i;
    CSpectrumFluxAxis FluxAxis(n, 1);
    std::string lsfType = "GaussianVariableWidth";

    TFloat64List widthList(n, 0.5);
    std::shared_ptr<TLSFArguments> args = std::make_shared<TLSFGaussianVarWidthArgs>(widthList,spcAxis);
    std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

    TFloat64List widthList2(n+1, 0.5);
    std::shared_ptr<TLSFArguments> args2 = std::make_shared<TLSFGaussianVarWidthArgs>(widthList2,spcAxis);
    BOOST_CHECK_EXCEPTION( LSFFactory.Create(lsfType, args2);,
                        std::exception, correctMessage);

}

BOOST_AUTO_TEST_CASE(LSFArgsPolymorphism)
{
    TFloat64List v = {1., 2., 3.}, w ={10., 20., 30.};
    std::shared_ptr<TLSFArguments> args = std::make_shared<TLSFGaussianVarWidthArgs>(v,w);
    //    TLSFGaussianVarWidthArgs *pB = dynamic_cast<TLSFGaussianVarWidthArgs*>(pA);

    std::string lsfType = "GaussianVariableWidth";
    std::shared_ptr<CLSF> LSF = LSFFactory.Create(lsfType, args);

    for(Int32 i=0; i<v.size(); i++)
    {
        BOOST_CHECK(LSF->GetWidth(v[i]) == w[i]);
    }

    BOOST_CHECK(LSF->GetWidth(1.2) == 12);
    //BOOST_CHECK(LSF->GetWidth(3.2) == 30);
    //BOOST_CHECK(LSF->GetWidth(0.2) == 10);

}

}
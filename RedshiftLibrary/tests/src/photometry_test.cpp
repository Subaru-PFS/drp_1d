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
#include "RedshiftLibrary/photometry/photometricdata.h"
#include "RedshiftLibrary/photometry/photometricband.h"

#include <exception>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;


BOOST_AUTO_TEST_SUITE(Photometry)

BOOST_AUTO_TEST_CASE(photometricdata)
{
    //const std::string name[] = {"band1", "band2"};
    const TStringList name = {"band1", "band2"};
    const TFloat64List flux = {1e-15, 1e-14};
    const TFloat64List fluxErr = {1e-18, 2e-18};


    BOOST_CHECK_NO_THROW(CPhotometricData phdat(name,flux,fluxErr));

    const TStringList name1 = {"band1"};
    BOOST_CHECK_THROW(CPhotometricData phdat(name1,flux,fluxErr),std::runtime_error);
    
    CPhotometricData phdat(name,flux,fluxErr);
    BOOST_CHECK(phdat.GetFlux("band1") == 1e-15);

    BOOST_CHECK(phdat.GetFluxErr("band2") == 2e-18);

    BOOST_CHECK(phdat.GetFluxOverErr2("band1") == 1e6);

    BOOST_CHECK(phdat.GetNameList() == name);

}

BOOST_AUTO_TEST_CASE(photometricbands)
{
    const Float64 trans[] = {.2, .5, .99, .99, .4, .1};
    const Float64 lambda[]= { 2000., 2200., 2321., 2430., 2554., 2603.};
    const Int32 n = sizeof(trans)/sizeof(trans[0]);
    
    BOOST_CHECK_NO_THROW(CPhotometricBand(trans,n,lambda,n));
    BOOST_CHECK_THROW(CPhotometricBand(trans,n,lambda,n+1), std::runtime_error);

    CPhotometricBand band(trans,n,lambda,n);
    BOOST_CHECK(band.m_transmission==TFloat64List(trans,trans+n));
    BOOST_CHECK(band.m_lambda==TFloat64List(lambda, lambda+n));

}

BOOST_AUTO_TEST_CASE(photometricbandcatalog)
{
    const Float64 trans1[] = {.2, .5, .99, .99, .4, .1};
    const Float64 lambda1[]= { 2000., 2200., 2321., 2430., 2554., 2603.};
    const Int32 n1 = sizeof(trans1)/sizeof(trans1[0]);

    const Float64 trans2[] = {.25, .55, .995, .995, .45, .15};
    const Float64 lambda2[]= { 4000., 4200., 4321., 4430., 4554., 4603.};
    const Int32 n2 = sizeof(trans2)/sizeof(trans2[0]);

    CPhotometricBand band1(trans1,n1,lambda1,n1);
    CPhotometricBand band2(trans2,n2,lambda2,n2);

    BOOST_CHECK_NO_THROW(CPhotBandCatalog  cat);

    CPhotBandCatalog  cat;
    TStringList names = {"band1", "band2"}; 
    BOOST_CHECK_NO_THROW(cat.Add(names[0], band1));
    BOOST_CHECK_NO_THROW(cat.Add(names[1], band2));
    BOOST_CHECK(cat.GetNameList() == names);
    BOOST_CHECK(cat.at(names[1]).m_transmission==TFloat64List(trans2,trans2+n2));
    BOOST_CHECK(cat.at(names[1]).m_lambda==TFloat64List(lambda2,lambda2+n2));

}

BOOST_AUTO_TEST_SUITE_END()

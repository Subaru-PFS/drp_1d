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
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/tests/test-tools.h"

#include <boost/test/unit_test.hpp>

using namespace NSEpic;
using namespace CPFTest;

BOOST_AUTO_TEST_SUITE(TemplateCatalog)

BOOST_AUTO_TEST_CASE(LoadCatalog)
{
    CTemplateCatalog catalog_w;
    CTemplateCatalog catalog_r;
    TStringList categories;

    boost::filesystem::path _path = boost::filesystem::unique_path("tst_%%%%%%%%%%");

    BOOST_REQUIRE(boost::filesystem::create_directories(_path));

    generate_template_catalog(catalog_r, 100, 3500., 12500.);

    BOOST_CHECK( catalog_r.GetTemplateCount( "galaxy" ) == 3 );
    BOOST_CHECK( catalog_r.GetTemplateCount( "emission" ) == 2 );
    BOOST_CHECK( catalog_r.GetTemplateCount( "qso" ) == 1 );
    BOOST_CHECK( catalog_r.GetTemplateCount( "star" ) == 2 );
    
    categories.push_back("galaxy");
    categories.push_back("star");
    TStringList expected={ "galaxy_test_template_2.txt",
                           "galaxy_test_template_1.txt",
                           "galaxy_test_template_0.txt",
                           "star_test_template_0.txt",
                           "star_test_template_1.txt"};
    bool found;
    const TTemplateRefList tplRef = catalog_r.GetTemplateList(categories);

    //BOOST_CHECK(expected.size() == tplRef.size());

    // look up expected template names in catalog
    for (UInt32 i=0; i<expected.size(); i++)
    {
        found = false;
        for (UInt32 j=0; j<tplRef.size(); j++)
        {
            if (tplRef[j]->GetName() == expected[i])
            {
                found = true;
                break;
            }
        }
        //BOOST_CHECK(found);
    }


    boost::filesystem::remove_all(_path);
}

BOOST_AUTO_TEST_SUITE_END()


#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/tests/test-tools.h>

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

    generate_template_catalog(catalog_w, 100, 3500., 12500.);

    BOOST_CHECK_NO_THROW(catalog_w.Save(_path.c_str(), true));

    BOOST_CHECK_THROW(catalog_r.Load( "/path/should/not/exist" ), std::runtime_error);

    BOOST_CHECK_THROW(catalog_r.Add( "/add/bogus/template/file", "galaxy" ), std::runtime_error);

    BOOST_CHECK_NO_THROW(catalog_r.Load( _path.c_str() ));

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
    TTemplateRefList tplRef = catalog_r.GetTemplate(categories);

    BOOST_CHECK(expected.size() == tplRef.size());

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
        BOOST_CHECK(found);
    }


    boost::filesystem::remove_all(_path);
}

BOOST_AUTO_TEST_SUITE_END()


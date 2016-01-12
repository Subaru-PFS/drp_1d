#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>


#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Template)


BOOST_AUTO_TEST_CASE(LoadCatalog)
{
    CTemplateCatalog catalog;

    Bool rValue = catalog.Load( "../test/data/templatecatalog/" );
    BOOST_CHECK( rValue == true );
    BOOST_CHECK( catalog.GetTemplateCount( "galaxy" ) == 2 );
    BOOST_CHECK( catalog.GetTemplateCount( "emission" ) == 1 );
    BOOST_CHECK( catalog.GetTemplateCount( "qso" ) == 1 );
    BOOST_CHECK( catalog.GetTemplateCount( "star" ) == 1 );

}


BOOST_AUTO_TEST_SUITE_END()


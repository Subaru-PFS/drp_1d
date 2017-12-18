#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>


#include <boost/test/unit_test.hpp>
#include "test-config.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Template)

BOOST_AUTO_TEST_CASE(LoadCatalog)
{
    CTemplateCatalog catalog;

    Bool rValue = catalog.Load( DATA_ROOT_DIR "templatecatalog/" );
    BOOST_CHECK( rValue == true );
    BOOST_CHECK( catalog.GetTemplateCount( "galaxy" ) == 2 );
    BOOST_CHECK( catalog.GetTemplateCount( "emission" ) == 1 );
    BOOST_CHECK( catalog.GetTemplateCount( "qso" ) == 1 );
    BOOST_CHECK( catalog.GetTemplateCount( "star" ) == 1 );

}

BOOST_AUTO_TEST_SUITE_END()


#include "template.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>

using namespace NSEpic;

using namespace NSEpicTest;

void CRedshiftTemplateTestCase::setUp()
{
}

void CRedshiftTemplateTestCase::tearDown()
{
}

void CRedshiftTemplateTestCase::LoadCatalog()
{
    CTemplateCatalog catalog;

    Bool rValue = catalog.Load( "../test/data/templatecatalog/" );
    CPPUNIT_ASSERT( rValue == true );
    CPPUNIT_ASSERT( catalog.GetTemplateCount( "galaxy" ) == 2 );
    CPPUNIT_ASSERT( catalog.GetTemplateCount( "emission" ) == 1 );
    CPPUNIT_ASSERT( catalog.GetTemplateCount( "qso" ) == 1 );
    CPPUNIT_ASSERT( catalog.GetTemplateCount( "star" ) == 1 );

}

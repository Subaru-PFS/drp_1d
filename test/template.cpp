#include "template.h"

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>

using namespace NSEpic;

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
    CPPUNIT_ASSERT( catalog.GetTemplateCount( CTemplate::nCategory_Galaxy ) == 2 );
    CPPUNIT_ASSERT( catalog.GetTemplateCount( CTemplate::nCategory_Emission ) == 1 );
    CPPUNIT_ASSERT( catalog.GetTemplateCount( CTemplate::nCategory_Qso ) == 1 );
    CPPUNIT_ASSERT( catalog.GetTemplateCount( CTemplate::nCategory_Star ) == 1 );

}

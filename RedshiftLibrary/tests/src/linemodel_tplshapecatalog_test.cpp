#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalogsTplShape.h>

#include <boost/test/unit_test.hpp>
#include "test-config.h"


using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(LinemodelTplshapecatalog)

BOOST_AUTO_TEST_CASE( LoadCatalogAndVelocities )
{
    std::string calibrationPath = DATA_ROOT_DIR "LinemodelTplshapeCtlgTestCase/calibrationPath1";
    std::string relPath = "linecatalogs_tplshape_ExtendedTemplatesMarch2016_B13B_mod20170110_3";
    CRayCatalogsTplShape* catalogTplShape = new CRayCatalogsTplShape();
    bool retVal = catalogTplShape->Init(calibrationPath, relPath);
    BOOST_CHECK( retVal == true );

    for(Int32 k=0; k<catalogTplShape->GetCatalogsCount(); k++)
    {
        std::string name = catalogTplShape->GetCatalogName(k);

        Float64 elv=-1;
        Float64 alv=-1;
        bool ret = catalogTplShape->GetCatalogVelocities(k, elv, alv);
        BOOST_CHECK( retVal == true );

        bool noCorrespondence = true;
        std::size_t foundstra;
        foundstra = name.find("vvds-reddestdataExtensionData");
        if (!(foundstra==std::string::npos)){
            BOOST_CHECK_MESSAGE( 780.0 == elv, "Tplshape-ctlg-velocities: check vvds-reddest EL failed" );
            BOOST_CHECK_MESSAGE( 480.0 == alv, "Tplshape-ctlg-velocities: check vvds-reddest ABS failed" );
            noCorrespondence = false;
        }
        foundstra = name.find("Rebinned-NEW-E-extendeddataExtensionData");
        if (!(foundstra==std::string::npos)){
            BOOST_CHECK_MESSAGE( 580.0 == elv, "Tplshape-ctlg-velocities: check Rebinned-NEW-E EL failed" );
            BOOST_CHECK_MESSAGE( 520.0 == alv, "Tplshape-ctlg-velocities: check Rebinned-NEW-E ABS failed" );
            noCorrespondence = false;
        }
        foundstra = name.find("COMBINE-ave-BX-highblue-AND-Scd");
        if (!(foundstra==std::string::npos)){
            BOOST_CHECK_MESSAGE( 360.0 == elv, "Tplshape-ctlg-velocities: check COMBINE-ave-BX-highblue-AND-Scd EL failed" );
            BOOST_CHECK_MESSAGE( 620.0 == alv, "Tplshape-ctlg-velocities: check COMBINE-ave-BX-highblue-AND-Scd ABS failed" );
            noCorrespondence = false;
        }

        BOOST_CHECK_MESSAGE( !noCorrespondence, "Tplshape-ctlg-velocities: no correspondence found for: "<< name.c_str() );

    }
    //BOOST_CHECK_MESSAGE( ratioIN == ratioMOD, "RatioRange: 3rd spectrum test failed = rule should not have been applied in this case" );
}

BOOST_AUTO_TEST_SUITE_END()

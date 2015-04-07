#include <epic/core/test/helper.h>

int main( int argc, char **argv )
{
    NSEpic::CTestHelper testHelper( argc, argv );

    if( testHelper.Configure( argc, argv ) )
        return testHelper.Run( "../test/reports/report.xml" );

    return false;
}


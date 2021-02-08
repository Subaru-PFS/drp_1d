#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/templatefittinglog.h>

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "test-config.h"
#include <math.h>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(TFLOG)

BOOST_AUTO_TEST_CASE(FindTplSpectralIndex)
{
    COperatorTemplateFittingLog tflog;
    const TAxisSampleList  spcLambda={1., 2., 3., 4., 5.};
    const TAxisSampleList  tplLambda={1.2, 2.2, 3.3, 4.3, 5.3, 6.3, 7.4, 8.4, 9.5, 10.5};
    TFloat64Range range(0,5);
    Float64 step = 0.01;
    TInt32Range result, ref(0, 7);


    /*result = tflog.FindTplSpectralIndex( range, step, spcLambda, tplLambda);
    BOOST_CHECK(result.GetBegin()== ref.GetBegin());
    BOOST_CHECK(result.GetEnd()== ref.GetEnd());*/
}

BOOST_AUTO_TEST_SUITE_END()
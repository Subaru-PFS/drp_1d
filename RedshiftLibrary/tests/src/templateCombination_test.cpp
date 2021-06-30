#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/operator/tplcombination.h"

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include "test-config.h"
#include <math.h>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(TemplateCombination)

BOOST_AUTO_TEST_CASE(inverseMatrix)
{
/*
    UInt32 dim = 3;
    cov = gsl_matrix_alloc (dim, dim);
    TFloatList _cov(1.0, 0.6, 0.0,
                    0.0, 1.5, 1.0,
                    0.0, 1.0, 1.0 );
    Int32 count = 0;
    for(Int32 i = 0; i < dim; i++)
    {
        for(Int32 j = 0; j < dim; j++)
        {
            gsl_matrix_set (cov, i, j, _cov[count]);
            count++;
        }
    }

    //COperatorTplcombination::InverseMatrix(cov);
    TFloatList inv_cov(
            1,  -1.2,   1.2,
            0,   2.0,  -2.0, 
            0,  -2.0,   3.0);
    count = 0;
    for (Int32 i = 0; i < dim; ++i)
    {
        for (Int32 j = 0; j < dim; ++j)
        {
            Float64 v = gsl_matrix_get(&inv.matrix,i,j));
            //BOOST_CHECK(v == inv_cov[count]);
            count++;
        }
    }
   */      

}

BOOST_AUTO_TEST_SUITE_END()
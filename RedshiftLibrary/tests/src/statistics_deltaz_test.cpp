#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/chisquareresult.h>

#include <RedshiftLibrary/statistics/deltaz.h>

#include <fstream>
#include <boost/math/special_functions.hpp>

//#include <boost/format.hpp>
//#include <boost/filesystem/fstream.hpp>
//#include <boost/filesystem.hpp>
//#include <boost/filesystem/fstream.hpp>
//#include <boost/algorithm/string/predicate.hpp>

#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Statistics_deltaz)


void DeltazTestCompute( const char* chi2Path, const Float64 redshift, const TFloat64Range range )
{
    std::ifstream input( chi2Path );
    BOOST_CHECK( input.is_open() );

    CChisquareResult chi2Result;
    chi2Result.Load( input );
    BOOST_CHECK_MESSAGE( chi2Result.ChiSquare.size()>0, "Loaded Chisquare result is empty." );

    CDeltaz deltaz;
    Float64 dz=-1.;
    Int32 ret = deltaz.Compute(chi2Result.ChiSquare, chi2Result.Redshifts, redshift, range, dz);
    BOOST_CHECK_MESSAGE( ret==0, "deltaz process returned error" );

    BOOST_CHECK_MESSAGE( dz>0.5e-4 && dz<3e-4, "Deltaz is nul or negative" );
}

BOOST_AUTO_TEST_CASE(Deltaz)
{
    //
    Float64 center_redshift = 2.6238;
    Float64 zRangeHalf = 0.005;
    TFloat64Range redshiftRange = TFloat64Range( center_redshift-zRangeHalf, center_redshift+zRangeHalf );

    DeltazTestCompute( "./RedshiftLibrary/tests/src/data/DeltazTestCase/simulm201605_tplshapeconttplfit_linemodelsolve.linemodel.csv", center_redshift, redshiftRange);

}


BOOST_AUTO_TEST_SUITE_END()

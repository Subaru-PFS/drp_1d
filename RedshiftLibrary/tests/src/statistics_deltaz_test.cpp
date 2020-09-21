#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/processflow/resultstore.h>
#include <RedshiftLibrary/statistics/deltaz.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <istream>
#include <iostream>
#include <fstream>
#include <boost/math/special_functions.hpp>
#include <boost/filesystem.hpp>

#include <boost/test/unit_test.hpp>
using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Statistics_deltaz)

string chi2sample = "#Redshifts	ChiSquare	Overlap\n"
    "0	                                4.98333524134186914307065308094025e+04\n"
    "0.00010000000000000000479217360239	4.98333424134186934679746627807617e+04\n"
    "0.00020000000000000000958434720477	4.98333524134186955052427947521210e+04\n"
    "0.00030000000000000002792904796323	4.98333624134186975425109267234802e+04\n"
    "0.00040000000000000001916869440954	4.98333724134186995797790586948395e+04\n"
    "0.00050000000000000001040834085586	4.98333824134187016170471906661987e+04\n"
    "0.00060000000000000005585809592645	4.98333924134187036543153226375580e+04\n"
    "0.00069999999999999999288763374850	4.98334024134187056915834546089172e+04\n"
    "0.00080000000000000003833738881909	4.98334124134187077288515865802765e+04\n"
    "0.00090000000000000008378714388968	4.98334224134187097661197185516357e+04\n"
    "0.00100000000000000002081668171172	4.98334524134187118033878505229950e+04\n"
    "0.00110000000000000006626643678231	4.98334424134187138406559824943542e+04\n"
    "0.00120000000000000011171619185291	4.98334524134187158779241144657135e+04\n"
    "0.00130000000000000015716594692350	4.98334624134187179151922464370728e+04\n"
    "0.00139999999999999998577526749699	4.98334724134187199524603784084320e+04\n"
    "0.00150000000000000003122502256758	4.98334824134187219897285103797913e+04\n"
    "0.00160000000000000007667477763817	4.98334924134187240269966423511505e+04\n"
    "0.00170000000000000012212453270877	4.98335024134187260642647743225098e+04\n"
    "0.00180000000000000016757428777936	4.98335124134187281015329062938690e+04\n"
    "0.00189999999999999999618360835285	4.98335224134187301388010382652283e+04\n"
    "0.00200000000000000004163336342344	4.94605935174078622367233037948608e+04\n"
    "0.00210000000000000030392355299114	4.94606035174078642739914357662201e+04\n"
    "0.00220000000000000013253287356463	4.94606135174078663112595677375793e+04\n"
    "0.00229999999999999996114219413812	4.94606235174078683485276997089386e+04\n"
    "0.00240000000000000022343238370581	4.94606335174078703857958316802979e+04\n"
    "0.00250000000000000005204170427930	4.94606435174078724230639636516571e+04\n"
    "0.00260000000000000031433189384700	4.94606535174078744603320956230164e+04\n"
    "0.00270000000000000014294121442049	4.94606635174078764976002275943756e+04\n"
    "0.00279999999999999997155053499398	4.94606735174078785348683595657349e+04\n"
    "0.00290000000000000023384072456167	4.94606835174078805721364915370941e+04\n"
    "0.00300000000000000006245004513517	4.94273189573849595035426318645477e+04\n"
    "0.00310000000000000032474023470286	4.94273289573849615408107638359070e+04\n"
    "0.00320000000000000015334955527635	4.94273389573849635780788958072662e+04\n"
    "0.00329999999999999998195887584984	4.94273489573849656153470277786255e+04\n"
    "0.00340000000000000024424906541753	4.94273589573849676526151597499847e+04\n"
    "0.00350000000000000007285838599103	4.94273689573849696898832917213440e+04\n"
    "0.00360000000000000033514857555872	4.94273789573849717271514236927032e+04\n"
    "0.00370000000000000016375789613221	4.94373889573849737644195556640625e+04\n"
    "0.00379999999999999999236721670570	4.94573989573849758016876876354218e+04\n";

TFloat64Range redshiftRange = TFloat64Range( 0, 5);
Float64       redshiftStep = 1E-4;
TRedshiftList pdfz = redshiftRange.SpreadOverLog( redshiftStep );

void DeltazTestCompute( const string& sample, const Float64 redshift, const TFloat64Range range )
{
    istringstream input( sample );
    boost::filesystem::path temp = boost::filesystem::unique_path();
    ofstream output(temp.native());

    CChisquareResult chi2Result;
    chi2Result.Load( input );
    BOOST_CHECK_MESSAGE( chi2Result.ChiSquare.size()>0, "Loaded Chisquare result is empty." );

    COperatorResultStore result_store;
    CParameterStore param_store;
    CDataStore dummy(result_store, param_store);

    BOOST_CHECK_NO_THROW( chi2Result.Save(dummy, output) );

    boost::filesystem::remove(temp);

    CDeltaz deltaz;
    Float64 dz=-1.;
    Int32 iz, izmin, izmax;
    Int32 ret = deltaz.GetRangeIndices(chi2Result.Redshifts, redshift, range, iz, izmin, izmax );
    ret = deltaz.Compute(chi2Result.ChiSquare, chi2Result.Redshifts, iz, izmin, izmax, dz);
    BOOST_CHECK_MESSAGE( ret==0, "deltaz process returned error" );
    BOOST_CHECK_MESSAGE( dz>0, "Deltaz is nul or negative : " << dz );
}

BOOST_AUTO_TEST_CASE(Deltaz)
{
    //
    Float64 center_redshift = 0.001;
    Float64 zRangeHalf = 0.0001;
    TFloat64Range redshiftRange = TFloat64Range( center_redshift-zRangeHalf, center_redshift+zRangeHalf );

    DeltazTestCompute( chi2sample, center_redshift, redshiftRange);

}

BOOST_AUTO_TEST_CASE(Deltazbordermin)
{
    //
    Float64 center_redshift = 0;
    Float64 zRangeHalf = 0.0001;
    TFloat64Range redshiftRange = TFloat64Range( center_redshift-zRangeHalf, center_redshift+zRangeHalf );

    DeltazTestCompute( chi2sample, center_redshift, redshiftRange);

}
BOOST_AUTO_TEST_CASE(Deltazbordermax)
{
    //
    Float64 center_redshift = 0.0038;
    Float64 zRangeHalf = 0.0001;
    TFloat64Range redshiftRange = TFloat64Range( center_redshift-zRangeHalf, center_redshift+zRangeHalf );

    DeltazTestCompute( chi2sample, center_redshift, redshiftRange);

}
//both redshifts belong to overlapping range
BOOST_AUTO_TEST_CASE(Deltaz_overlapping1)
{
    TRedshiftList center_redshifts = {1.0, 1.5};
    TRedshiftList deltaz = { 0.5/3, 0.5/3};
    TFloat64RangeList ranges;
    TFloat64RangeList correct_ranges = {{0.5, 1.24990}, {1.25, 2}};

    CPdfCandidateszResult res;
    res.Init(center_redshifts, deltaz);
    res.SetIntegrationWindows(pdfz, ranges);

    BOOST_CHECK_CLOSE(ranges[0].GetEnd(), correct_ranges[0].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetEnd(), correct_ranges[1].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[0].GetBegin(), correct_ranges[0].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetBegin(), correct_ranges[1].GetBegin(), 1E-4);
}
//no overlapping
BOOST_AUTO_TEST_CASE(Deltaz_nooverlapping)
{
    TRedshiftList center_redshifts = {1.0, 4.0};
    TRedshiftList deltaz = { 0.5/3, 0.5/3};
    TFloat64RangeList ranges;
    TFloat64RangeList correct_ranges = {{0.5, 1.5}, {3.5, 4.5}};

    CPdfCandidateszResult res; 
    res.Init(center_redshifts, deltaz);
    res.SetIntegrationWindows(pdfz, ranges);

    BOOST_CHECK_CLOSE(ranges[0].GetEnd(), correct_ranges[0].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetEnd(), correct_ranges[1].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[0].GetBegin(), correct_ranges[0].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetBegin(), correct_ranges[1].GetBegin(), 1E-4);
}
//only the smallest redshift belongs to the overlapping region
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_3)
{
    TRedshiftList center_redshifts = {1.0, 4.0};
    TRedshiftList deltaz = { 1./3, 2.5/3};
    TFloat64RangeList ranges;
    TFloat64RangeList correct_ranges = {{0., 1.7499}, {1.75, pdfz.back()}};

    CPdfCandidateszResult res; 
    res.Init(center_redshifts, deltaz);
    res.SetIntegrationWindows(pdfz, ranges);

    BOOST_CHECK_CLOSE(ranges[0].GetEnd(), correct_ranges[0].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetEnd(), correct_ranges[1].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[0].GetBegin(), correct_ranges[0].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetBegin(), correct_ranges[1].GetBegin(), 1E-4);
}

//only the greatest redshift belong to overlapping region
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_4)
{
    TRedshiftList center_redshifts = {1.0, 4.0};
    TRedshiftList deltaz = {2.0/3, 2.5/3};
    TFloat64RangeList ranges;
    TFloat64RangeList correct_ranges = {{0., 2.2499}, {2.25, pdfz.back()}};

    CPdfCandidateszResult res;
    res.Init(center_redshifts, deltaz);
    res.SetIntegrationWindows(pdfz, ranges);

    BOOST_CHECK_CLOSE(ranges[0].GetEnd(), correct_ranges[0].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetEnd(), correct_ranges[1].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[0].GetBegin(), correct_ranges[0].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges[1].GetBegin(), correct_ranges[1].GetBegin(), 1E-4);
}
BOOST_AUTO_TEST_SUITE_END()

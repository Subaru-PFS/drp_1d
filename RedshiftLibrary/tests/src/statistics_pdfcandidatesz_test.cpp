#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/statistics/pdfcandidatesz.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"

#include <istream>
#include <iostream>
#include <fstream>
#include <boost/math/special_functions.hpp>
#include <boost/filesystem.hpp>

#include <boost/test/unit_test.hpp>
using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Statistics_pdfcandidatesz)

TFloat64Range redshiftRange = TFloat64Range( 0, 5);
Float64       redshiftStep = 1E-4;
TRedshiftList pdfz = redshiftRange.SpreadOverLogZplusOne( redshiftStep );



//both redshifts belong to overlapping range
BOOST_AUTO_TEST_CASE(Deltaz_overlapping1)
{
    TRedshiftList center_redshifts = {1.0, 1.5};
    TRedshiftList deltaz = { 0.5/3, 0.5/3};
    TCandidateZRangebyID ranges;
    TCandidateZRangebyID correct_ranges = { {"EXT0", {0.5, 1.24990}}, {"EXT1", {1.25, 2}} };

    TCandidateZbyID zcandidates;
    zcandidates["EXT0"].Redshift = center_redshifts[0];
    zcandidates["EXT0"].Deltaz = deltaz[0];
    zcandidates["EXT1"].Redshift = center_redshifts[1];
    zcandidates["EXT1"].Deltaz = deltaz[1];
    CPdfCandidatesZ zcand_op =  CPdfCandidatesZ(zcandidates);
    zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

    BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(), correct_ranges["EXT0"].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(), correct_ranges["EXT1"].GetBegin(), 1E-4);
}

//no overlapping
BOOST_AUTO_TEST_CASE(Deltaz_nooverlapping)
{
    TRedshiftList center_redshifts = {1.0, 4.0};
    TRedshiftList deltaz = { 0.5/3, 0.5/3};
    TCandidateZRangebyID ranges;
    TCandidateZRangebyID correct_ranges = { {"EXT0",{0.5, 1.5}}, {"EXT1", {3.5, 4.5}} };

    TCandidateZbyID zcandidates;
    zcandidates["EXT0"].Redshift = center_redshifts[0];
    zcandidates["EXT0"].Deltaz = deltaz[0];
    zcandidates["EXT1"].Redshift = center_redshifts[1];
    zcandidates["EXT1"].Deltaz = deltaz[1];
    CPdfCandidatesZ zcand_op =  CPdfCandidatesZ(zcandidates);
    zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

    BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(), correct_ranges["EXT0"].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(), correct_ranges["EXT1"].GetBegin(), 1E-4);
}
//only the smallest redshift belongs to the overlapping region
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_3)
{
    TRedshiftList center_redshifts = {1.0, 4.0};
    TRedshiftList deltaz = { 1./3, 2.5/3};
    TCandidateZRangebyID ranges;
    TCandidateZRangebyID correct_ranges = { {"EXT0", {0., 1.7499}} , {"EXT1", {1.75, pdfz.back()}} };

    TCandidateZbyID zcandidates;
    zcandidates["EXT0"].Redshift = center_redshifts[0];
    zcandidates["EXT0"].Deltaz = deltaz[0];
    zcandidates["EXT1"].Redshift = center_redshifts[1];
    zcandidates["EXT1"].Deltaz = deltaz[1];
    CPdfCandidatesZ zcand_op =  CPdfCandidatesZ(zcandidates);
    zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

    BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(), correct_ranges["EXT0"].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(), correct_ranges["EXT1"].GetBegin(), 1E-4);
}

//only the greatest redshift belong to overlapping region
BOOST_AUTO_TEST_CASE(Deltaz_overlapping_4)
{
    TRedshiftList center_redshifts = {1.0, 4.0};
    TRedshiftList deltaz = {2.0/3, 2.5/3};
    TCandidateZRangebyID ranges;
    TCandidateZRangebyID correct_ranges = { {"EXT0", {0., 2.2499}}, {"EXT1", {2.25, pdfz.back()}} };

    TCandidateZbyID zcandidates;
    zcandidates["EXT0"].Redshift = center_redshifts[0];
    zcandidates["EXT0"].Deltaz = deltaz[0];
    zcandidates["EXT1"].Redshift = center_redshifts[1];
    zcandidates["EXT1"].Deltaz = deltaz[1];
    CPdfCandidatesZ zcand_op =  CPdfCandidatesZ(zcandidates);
    zcand_op.SetIntegrationWindows(TFloat64Range(pdfz), ranges);

    BOOST_CHECK_CLOSE(ranges["EXT0"].GetEnd(), correct_ranges["EXT0"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetEnd(), correct_ranges["EXT1"].GetEnd(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT0"].GetBegin(), correct_ranges["EXT0"].GetBegin(), 1E-4);
    BOOST_CHECK_CLOSE(ranges["EXT1"].GetBegin(), correct_ranges["EXT1"].GetBegin(), 1E-4);
}
/*

BOOST_AUTO_TEST_CASE(SortByValSumProbaInt)
{
    TRedshiftList center_redshifts = {1.0, 4.0, 3.0};

    TCandidateZbyID zcandidates;
    zcandidates["EXT0"].Redshift = center_redshifts[1];
    zcandidates["EXT0"].ValSumProba = 0.;

    zcandidates["EXT1"].Redshift = center_redshifts[0];
    zcandidates["EXT1"].ValSumProba = 0.;

    zcandidates["EXT2"].Redshift = center_redshifts[2];
    zcandidates["EXT2"].ValSumProba = 1.;

    CPdfCandidatesZ zcand_op =  CPdfCandidatesZ(zcandidates);

    TCandidateZbyRank ranked_candidates;
    zcand_op.SortByValSumProbaInt(ranked_candidates);
    std::cout<< ranked_candidates[0].second.ValSumProba<<","<<ranked_candidates[0].first<<"\n";
    std::cout<< ranked_candidates[1].second.ValSumProba<<","<<ranked_candidates[1].first<<"\n";
    std::cout<< ranked_candidates[2].second.ValSumProba<<","<<ranked_candidates[2].first<<"\n";
}


BOOST_AUTO_TEST_CASE(SortByValSumProbaInt2)
{
    TCandidateZbyID zcandidates;
    zcandidates["EXT5"].ValSumProba = 0.;
    zcandidates["EXT1"].ValSumProba = 0.;
    zcandidates["EXT2"].ValSumProba = 1.;
    std::vector<std::string> Ids;
    for (const auto & c : zcandidates){
        Ids.push_back(c.first); // keys = ids
        std::cout<< c.first <<"\n";
    }
    const TCandidateZbyID & c = zcandidates; 

    std::stable_sort(Ids.rbegin(), Ids.rend(),
        [&c](std::string Id1, std::string Id2) {return c.at(Id1).ValSumProba < c.at(Id2).ValSumProba;});
    
    std::cout<<"Descending order after stable_sort: \n";
    TCandidateZbyRank ranked_candidates;
    for (const auto & Id: Ids){
        std::cout<<Id<<" "<<zcandidates[Id].ValSumProba<<"\n";
        ranked_candidates.emplace_back(Id, zcandidates.at(Id));
    }
    for (const auto & cc : ranked_candidates){
        std::cout<<cc.first<<" "<<cc.second.ValSumProba<<"\n";
    }
}*/
BOOST_AUTO_TEST_SUITE_END()

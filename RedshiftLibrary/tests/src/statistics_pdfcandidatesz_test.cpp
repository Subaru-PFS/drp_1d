// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
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

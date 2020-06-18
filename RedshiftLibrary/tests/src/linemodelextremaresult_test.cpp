#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/linemodel/linemodelextremaresult.h>
#include <iostream>
#include <boost/test/unit_test.hpp>

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(linemodelextremaresult)


BOOST_AUTO_TEST_CASE(Extremum1)
{
  CLineModelExtremaResult m_firstpass_extremaResult;
  m_firstpass_extremaResult.Extrema = {1.0, 2.0, 3.0, 4.0, 5.0};
  TInt32List Rank_PDF = {0, 2, 3, 1};
  std::vector<std::string> ids =  {"FPE0", "FPE2", "FPE4", "FPE3"};
  TInt32List Rank_PDF_firstpass = {0, 2, 4, 3, 1};
  Int32 ret = m_firstpass_extremaResult.FixRanksUsingSortedIDs( Rank_PDF, ids); 
  BOOST_CHECK(ret==0);
  BOOST_CHECK_EQUAL_COLLECTIONS(Rank_PDF.begin(), Rank_PDF.begin()+5, 
                              Rank_PDF_firstpass.begin(), Rank_PDF_firstpass.begin()+5);


}

BOOST_AUTO_TEST_SUITE_END()

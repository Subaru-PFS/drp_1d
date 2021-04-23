#include <RedshiftLibrary/linemodel/zweimodelresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <RedshiftLibrary/spectrum/spectrum.h>

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CZweiModelResult::CZweiModelResult()
{
  
}

CZweiModelResult::CZweiModelResult(std::vector<Float64> redshifts_s1, std::vector<Float64> redshifts_s2, std::vector<std::vector<Float64>> combined_merits)
{
    m_redshifts_s1 = redshifts_s1;
    m_redshifts_s2 = redshifts_s2;
    m_combined_merits = combined_merits;
}

/**
 * \brief Empty destructor.
 **/
CZweiModelResult::~CZweiModelResult()
{

}


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

void CZweiModelResult::Save(std::ostream& stream ) const
{
    // first save the s1 redshift values on 1 line
    stream << "#redhifts s1" << std::endl;;
    for(Int32 k=0; k<m_redshifts_s1.size(); k++)
    {
        stream << m_redshifts_s1[k] << "\t";
    }
    stream << std::endl;

    // then save the s2 redshift values on 1 line
    stream << "#redhifts s2" << std::endl;;
    for(Int32 k=0; k<m_redshifts_s2.size(); k++)
    {
        stream << m_redshifts_s2[k] << "\t";
    }
    stream << std::endl;

    // finally save the merits matrix (rows=s2 index, cols = s1 index)
    stream << "#merits matrix (rows=s2 index, cols = s1 index)" << std::endl;
    for(Int32 k2=0; k2<m_redshifts_s2.size(); k2++)
    {

        for(Int32 k1=0; k1<m_redshifts_s1.size(); k1++)
        {
            stream << m_combined_merits[k2][k1] << "\t";
        }
        stream << std::endl;
    }
}

/**
 * \brief Empty method.
 **/
void CZweiModelResult::SaveLine( std::ostream& stream ) const
{
    stream << "testSaveLine";
}


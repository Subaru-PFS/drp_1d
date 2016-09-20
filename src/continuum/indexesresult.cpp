
#include <epic/redshift/continuum/indexesresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <epic/redshift/spectrum/spectrum.h>

using namespace NSEpic;


/**
 * \brief Empty constructor.
 **/
CContinuumIndexesResult::CContinuumIndexesResult()
{
}

/**
 * \brief Empty destructor.
 **/
CContinuumIndexesResult::~CContinuumIndexesResult()
{
}


/**
 * \brief Set the indexes values
 **/
void CContinuumIndexesResult::SetValues(Float64 stdSpc, Float64 std_continuum)
{
    m_StdSpectrum = stdSpc;
    m_StdContinuum = std_continuum;
}

/**
 * \brief Prints the result in the argument store, using the argument stream as output.
 **/
Void CContinuumIndexesResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    stream << "std_spc" << "\t" << "std_continuum" << std::endl;
    stream << m_StdSpectrum << "\t" << m_StdContinuum << std::endl;
}

/**
 * \brief Empty method.
 **/
Void CContinuumIndexesResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}

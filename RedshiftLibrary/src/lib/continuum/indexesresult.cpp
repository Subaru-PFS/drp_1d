
#include <RedshiftLibrary/continuum/indexesresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <RedshiftLibrary/spectrum/spectrum.h>

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


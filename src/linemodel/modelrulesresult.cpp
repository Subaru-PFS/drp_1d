#include <epic/redshift/linemodel/modelrulesresult.h>

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
CModelRulesResult::CModelRulesResult(TStringList logStrings)
{
    LogStrings = logStrings;
}

/**
 * \brief Empty constructor.
 **/
CModelRulesResult::CModelRulesResult()
{
}

/**
 * \brief Empty destructor.
 **/
CModelRulesResult::~CModelRulesResult()
{
}

/**
 * \brief Prints the rules of the Linemodel in the argument store, using the argument stream as output.
 **/
Void CModelRulesResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    for(Int32 k=0; k<LogStrings.size(); k++)
    {
        stream << k << "\t" << LogStrings[k].c_str() << std::endl;
    }
}

/**
 * \brief Empty method.
 **/
Void CModelRulesResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}

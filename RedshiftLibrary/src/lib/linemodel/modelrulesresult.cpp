#include <RedshiftLibrary/linemodel/modelrulesresult.h>

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
void CModelRulesResult::Save(std::ostream& stream ) const
{
    for(Int32 k=0; k<LogStrings.size(); k++)
    {
        stream << k << "\t" << LogStrings[k].c_str() << std::endl;
    }
}

/**
 * \brief Empty method.
 **/
void CModelRulesResult::SaveLine(std::ostream& stream ) const
{

}

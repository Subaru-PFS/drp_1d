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


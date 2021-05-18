#ifndef _REDSHIFT_METHOD_TPLCOMBINATIONSOLVERESULT_
#define _REDSHIFT_METHOD_TPLCOMBINATIONSOLVERESULT_

#include <RedshiftLibrary/method/solveresult.h>
#include <RedshiftLibrary/operator/extremaresult.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <memory>
#include <vector>
#include <unordered_map>
#include <cmath>

namespace NSEpic
{


/**
 * \ingroup Redshift
 */
class CTplCombinationSolveResult : public CPdfSolveResult
{

public:
    CTplCombinationSolveResult( Float64 merit, Float64 redshift,
                                const std::string & scope, 
                                const std::string & opt_pdfcombination,
                                Float64 evidence );

  ~CTplCombinationSolveResult();
private:

    const std::string m_scope;

};

}

#endif

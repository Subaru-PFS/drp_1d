#include "RedshiftLibrary/method/tplcombinationsolveresult.h"

#include <stdio.h>
#include <float.h>
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/method/solveresult.h"
#include <memory>

using namespace NSEpic;
CTplCombinationSolveResult::CTplCombinationSolveResult(const std::string & scope, 
                                                         const TCandidateZ & BestExtremumResult,
                                                         const std::string & opt_pdfcombination,
                                                         Float64 evidence ):
    CPdfSolveResult( BestExtremumResult, opt_pdfcombination, evidence),
    m_scope(scope)
{}

CTplCombinationSolveResult::~CTplCombinationSolveResult(){}

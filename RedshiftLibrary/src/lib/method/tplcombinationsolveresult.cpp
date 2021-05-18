#include <RedshiftLibrary/method/tplcombinationsolveresult.h>

#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/common/exception.h>
#include <RedshiftLibrary/common/formatter.h>
#include <RedshiftLibrary/method/solveresult.h>
#include <memory>

using namespace NSEpic;
CTplCombinationSolveResult::CTplCombinationSolveResult( Float64 merit, Float64 redshift,
                                                        const std::string & scope, 
                                                        const std::string & opt_pdfcombination,
                                                        Float64 evidence ):
    CPdfSolveResult( merit, redshift, opt_pdfcombination, evidence),
    m_scope(scope)
{}

CTplCombinationSolveResult::~CTplCombinationSolveResult(){}
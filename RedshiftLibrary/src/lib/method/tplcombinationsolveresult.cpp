#include <RedshiftLibrary/method/tplcombinationsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/extremum/extremum.h>

using namespace NSEpic;

CTplcombinationSolveResult::CTplcombinationSolveResult()
{
    m_type = nType_raw;
    m_scope = "tplcombination";
    m_name = "Tplcombination";
}

CTplcombinationSolveResult::~CTplcombinationSolveResult()
{

}


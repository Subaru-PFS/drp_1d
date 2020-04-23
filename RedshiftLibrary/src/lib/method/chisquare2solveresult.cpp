#include <RedshiftLibrary/method/chisquare2solveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>
#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/extremum/extremum.h>

using namespace NSEpic;

CChisquare2SolveResult::CChisquare2SolveResult()
{
    m_type = nType_raw;
    m_scope = "chisquare2";
    m_name = "Chisquare2";
}

CChisquare2SolveResult::~CChisquare2SolveResult()
{

}


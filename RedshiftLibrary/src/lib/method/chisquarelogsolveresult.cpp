#include <RedshiftLibrary/method/chisquarelogsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/extremum/extremum.h>


using namespace NSEpic;

CChisquareLogSolveResult::CChisquareLogSolveResult()
{
    m_type = nType_raw;
    m_scope = "chisquarelog";
    m_name = "ChisquareLog";
}

CChisquareLogSolveResult::~CChisquareLogSolveResult()
{

}


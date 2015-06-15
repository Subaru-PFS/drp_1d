#include <epic/redshift/operator/blindsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CBlindSolveResult )

CBlindSolveResult::CBlindSolveResult()
{

}

CBlindSolveResult::~CBlindSolveResult()
{

}

Void CBlindSolveResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestFitResult( store, redshift, merit, tplName );

    stream <<  "#Spectrum\tRedshifts\tMerit\tTemplate"<< std::endl;

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << std::endl;
}

Bool CBlindSolveResult::GetBestFitResult( const COperatorResultStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const
{
    TOperatorResultMap correlationResults = store.GetPerTemplateResult("blindsolve.correlation");
    TOperatorResultMap meritResults = store.GetPerTemplateResult("blindsolve.merit");


    Int32 maxIndex = 0;
    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        const CChisquareResult* meritResult = (const CChisquareResult*)(const COperatorResult*)(*it).second;
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->ChiSquare[i] < tmpMerit && meritResult->Status[i] == COperator::nStatus_OK )
            {
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpTplName = (*it).first;
            }
        }
    }


    if( tmpMerit < DBL_MAX )
    {
        redshift = tmpRedshift;
        merit = tmpMerit;
        tplName = tmpTplName;
        return true;
    }

    return false;

}

Bool CBlindSolveResult::GetBestCorrelationPeakResult( const COperatorResultStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const
{
    TOperatorResultMap correlationResults = store.GetPerTemplateResult("blindsolve.correlation");
    TOperatorResultMap meritResults = store.GetPerTemplateResult("blindsolve.merit");


    Int32 maxIndex = 0;
    Float64 tmpCorr = DBL_MIN ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName;

    for( TOperatorResultMap::const_iterator it = correlationResults.begin(); it != correlationResults.end(); it++ )
    {
        const CCorrelationResult* corrResult = (const CCorrelationResult*)(const COperatorResult*)(*it).second;
        for( Int32 i=0; i<corrResult->Correlation.size(); i++ )
        {
            if( corrResult->Correlation[i] > tmpCorr && corrResult->Status[i] == COperator::nStatus_OK )
            {
                tmpCorr = corrResult->Correlation[i];
                tmpRedshift = corrResult->Redshifts[i];
                tmpTplName = (*it).first;
            }
        }
    }


    if( tmpCorr > DBL_MIN )
    {
        redshift = tmpRedshift;
        merit = tmpCorr;
        tplName = tmpTplName;
        return true;
    }

    return false;

}

#include <epic/redshift/method/fullsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CFullSolveResult )

CFullSolveResult::CFullSolveResult()
{

}

CFullSolveResult::~CFullSolveResult()
{

}

Void CFullSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestCorrelationResult( store, redshift, merit, tplName );

    stream <<  "#Spectrum\tRedshifts\tMerit\tTemplate"<< std::endl;

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << std::endl;
}

Void CFullSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestCorrelationResult( store, redshift, merit, tplName );

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "FullSolve" << std::endl;
}


Bool CFullSolveResult::GetBestCorrelationResult( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const
{
    std::string scope = store.GetScope( this ) + "fullsolve.correlation";
    TOperatorResultMap correlationResults = store.GetPerTemplateResult(scope.c_str());


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

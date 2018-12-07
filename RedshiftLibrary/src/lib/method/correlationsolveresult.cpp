#include <RedshiftLibrary/method/correlationsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;


CCorrelationSolveResult::CCorrelationSolveResult()
{

}

CCorrelationSolveResult::~CCorrelationSolveResult()
{

}

void CCorrelationSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
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


void CCorrelationSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestCorrelationResult( store, redshift, merit, tplName );

    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "CorrelationSolve" << "\t"
                << "-1" << "\t" //deltaz
                << "-1" << std::endl; //reliability label
}


Bool CCorrelationSolveResult::GetBestCorrelationResult( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const
{
    std::string scope = store.GetScope( *this ) + "correlationsolve.correlation";
    TOperatorResultMap correlationResults = store.GetPerTemplateResult(scope.c_str());


    Float64 tmpCorr = DBL_MIN ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName;

    for( TOperatorResultMap::const_iterator it = correlationResults.begin(); it != correlationResults.end(); it++ )
    {
        auto corrResult = std::dynamic_pointer_cast<const CCorrelationResult>( (*it).second );
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

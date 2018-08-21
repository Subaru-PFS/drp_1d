#include <RedshiftLibrary/method/blindsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

CBlindSolveResult::CBlindSolveResult()
{

}

CBlindSolveResult::~CBlindSolveResult()
{

}

void CBlindSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
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

void CBlindSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestFitResult( store, redshift, merit, tplName );
    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "BlindSolve" << "\t"
                << "-1" << "\t" //deltaz
                << "-1" << std::endl; //reliability label
}

Bool CBlindSolveResult::GetBestFitResult( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const
{
    std::string scope_corr = store.GetScope( *this ) + "blindsolve.correlation";
    TOperatorResultMap correlationResults = store.GetPerTemplateResult( scope_corr.c_str() );
    //TOperatorResultMap correlationResults = store.GetPerTemplateResult("blindsolve.correlation");

    std::string scope_merit = store.GetScope( *this ) + "blindsolve.merit";
    TOperatorResultMap meritResults = store.GetPerTemplateResult( scope_merit.c_str() );


    Int32 maxIndex = 0;
    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast< const CChisquareResult>( (*it).second );
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

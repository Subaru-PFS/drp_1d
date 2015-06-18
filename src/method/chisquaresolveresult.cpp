#include <epic/redshift/method/chisquaresolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CChisquareSolveResult )

CChisquareSolveResult::CChisquareSolveResult()
{

}

CChisquareSolveResult::~CChisquareSolveResult()
{

}

Void CChisquareSolveResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit, tplName );

    stream <<  "#Spectrum\tRedshifts\tMerit\tTemplate"<< std::endl;

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << std::endl;
}

Void CChisquareSolveResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit, tplName );

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "ChisquareSolve" << std::endl;
}

Bool CChisquareSolveResult::GetBestRedshift( const COperatorResultStore& store, Float64& redshift, Float64& merit, std::string& tplName ) const
{

    std::string scope = store.GetScope( this ) + "chisquaresolve.chisquare";
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());


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



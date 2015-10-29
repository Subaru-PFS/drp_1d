#include <epic/redshift/method/linemodelsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/linemodelresult.h>
#include <stdio.h>
#include <float.h>
#include <math.h>


using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CLineModelSolveResult )

CLineModelSolveResult::CLineModelSolveResult()
{

}

CLineModelSolveResult::~CLineModelSolveResult()
{

}

Void CLineModelSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit );

    stream <<  "#Redshifts\tMerit\tTemplate"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << std::endl;


}


Void CLineModelSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit );

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "LineModelSolve" << std::endl;
}

Bool CLineModelSolveResult::GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit ) const
{

    std::string scope = store.GetScope( this ) + "linemodelsolve.linemodel";
    const CLineModelResult* results = (CLineModelResult*)store.GetGlobalResult(scope.c_str());


    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;

    if(results){
        for( Int32 i=0; i<results->ChiSquare.size(); i++ )
        {
            if( results->ChiSquare[i] < tmpMerit )
            {
                tmpMerit = results->ChiSquare[i];
                tmpRedshift = results->Redshifts[i];
            }
        }

    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;

}

Bool CLineModelSolveResult::GetBestRedshiftLogArea( const CDataStore& store, Float64& redshift, Float64& merit ) const
{

    std::string scope = store.GetScope( this ) + "linemodelsolve.linemodel";
    CLineModelResult* results = (CLineModelResult*)store.GetGlobalResult(scope.c_str());

    Float64 tmpMerit = -DBL_MAX ;
    Float64 tmpRedshift = 0.0;

    if(results){
        for( Int32 i=0; i<results->LogArea.size(); i++ )
        {
            if( results->LogArea[i] > tmpMerit )
            {
                tmpMerit = results->LogArea[i];
                tmpRedshift = results->LogAreaCorrectedExtrema[i];
            }
        }

    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
}




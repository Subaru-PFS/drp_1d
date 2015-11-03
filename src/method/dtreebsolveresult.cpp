#include <epic/redshift/method/dtreebsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <epic/redshift/operator/linemodelresult.h>

#include <float.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CDTreeBSolveResult )

CDTreeBSolveResult::CDTreeBSolveResult()
{

}

CDTreeBSolveResult::~CDTreeBSolveResult()
{

}

Void CDTreeBSolveResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( this ) + "dtreeBsolve.redshiftresult";
    const COperatorResult* Results = store.GetGlobalResult(scope.c_str());
    if(Results!=NULL){
        Results->Save(store, stream );
    }
}

Void CDTreeBSolveResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( this ) + "dtreeBsolve.redshiftresult";
    const COperatorResult* Results = store.GetGlobalResult(scope.c_str());
    if(Results!=NULL){
        Results->SaveLine(store, stream);
    }
}



Bool CDTreeBSolveResult::GetBestRedshift( const COperatorResultStore& store, Float64& redshift, Float64& merit ) const
{

    std::string scope = store.GetScope( this ) + "dtreeBsolve.chisquare2solve.chisquare";
    //std::string scope = store.GetScope( this ) + "dtreeBsolve.linemodel";
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

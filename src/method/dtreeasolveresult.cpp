#include <epic/redshift/method/dtreeasolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CDTreeASolveResult )

CDTreeASolveResult::CDTreeASolveResult()
{

}

CDTreeASolveResult::~CDTreeASolveResult()
{

}

Void CDTreeASolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( this ) + "dtreeAsolve.redshiftresult";
    const COperatorResult* Results = store.GetGlobalResult(scope.c_str());
    if(Results!=NULL){
        Results->Save(store, stream );
    }
}

Void CDTreeASolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( this ) + "dtreeAsolve.redshiftresult";
    const COperatorResult* Results = store.GetGlobalResult(scope.c_str());
    if(Results!=NULL){
        Results->SaveLine(store, stream);
    }
}




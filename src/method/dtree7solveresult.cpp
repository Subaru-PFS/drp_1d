#include <epic/redshift/method/dtree7solveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CDTree7SolveResult )

CDTree7SolveResult::CDTree7SolveResult()
{

}

CDTree7SolveResult::~CDTree7SolveResult()
{

}

Void CDTree7SolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( this ) + "dtree7solve.redshiftresult";
    const COperatorResult* Results = store.GetGlobalResult(scope.c_str());
    Results->Save(store, stream );
}

Void CDTree7SolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( this ) + "dtree7solve.redshiftresult";
    const COperatorResult* Results = store.GetGlobalResult(scope.c_str());
    Results->SaveLine(store, stream);
}





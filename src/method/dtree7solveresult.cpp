#include <epic/redshift/method/dtree7solveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

CDTree7SolveResult::CDTree7SolveResult()
{

}

CDTree7SolveResult::~CDTree7SolveResult()
{

}

Void CDTree7SolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( *this ) + "dtree7solve.redshiftresult";
    auto Results = store.GetGlobalResult(scope.c_str());
    Results.lock()->Save(store, stream );
}

Void CDTree7SolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( *this ) + "dtree7solve.redshiftresult";
    auto Results = store.GetGlobalResult(scope.c_str());
    Results.lock()->SaveLine(store, stream);
}





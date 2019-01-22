#include <RedshiftLibrary/method/dtree7solveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;

CDTree7SolveResult::CDTree7SolveResult()
{

}

CDTree7SolveResult::~CDTree7SolveResult()
{

}

void CDTree7SolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( *this ) + "dtree7solve.redshiftresult";
    auto Results = store.GetGlobalResult(scope.c_str());
    Results.lock()->Save(store, stream );
}

void CDTree7SolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( *this ) + "dtree7solve.redshiftresult";
    auto Results = store.GetGlobalResult(scope.c_str());
    Results.lock()->SaveLine(store, stream);
}





#include <RedshiftLibrary/method/dtreeasolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>

#include <float.h>

using namespace NSEpic;


CDTreeASolveResult::CDTreeASolveResult()
{

}

CDTreeASolveResult::~CDTreeASolveResult()
{

}

Void CDTreeASolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( *this ) + "dtreeAsolve.redshiftresult";
    auto Results = store.GetGlobalResult(scope.c_str());
    if(!Results.expired()){
        Results.lock()->Save(store, stream );
    }
}

Void CDTreeASolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    std::string scope = store.GetScope( *this ) + "dtreeAsolve.redshiftresult";
    auto Results = store.GetGlobalResult(scope.c_str());
    if(!Results.expired()){
        Results.lock()->SaveLine(store, stream);
    }
}




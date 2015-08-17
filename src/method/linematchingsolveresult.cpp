#include <epic/redshift/method/linematchingsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <float.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CLineMatchingSolveResult )

CLineMatchingSolveResult::CLineMatchingSolveResult()
{

}

CLineMatchingSolveResult::~CLineMatchingSolveResult()
{

}

Void CLineMatchingSolveResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;

    GetBestResult( store, redshift, merit );

    stream <<  "#Spectrum\tRedshifts\tMatchNum\t"<< std::endl;

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << std::endl;
}

Void CLineMatchingSolveResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;

    GetBestResult( store, redshift, merit );
    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << "LineMatchingSolve" << std::endl;
}

Bool CLineMatchingSolveResult::GetBestResult(const COperatorResultStore& store, Float64& redshift, Float64& merit) const
{
    std::string scope = store.GetScope( this ) + "linematchingsolve.raymatching";
    const CRayMatchingResult* Results = (CRayMatchingResult*)store.GetGlobalResult(scope.c_str());

    Int32 tmpMatchNum = -1;
    Float64 tmpRedshift = -1.0;

    if(Results){
        Int32 er = Results->GetBestRedshift(tmpRedshift, tmpMatchNum);
    }

    redshift = tmpRedshift;
    merit = tmpMatchNum;
    return true;
}

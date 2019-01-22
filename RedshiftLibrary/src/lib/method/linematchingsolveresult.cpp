#include <RedshiftLibrary/method/linematchingsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>

#include <float.h>

using namespace NSEpic;

CLineMatchingSolveResult::CLineMatchingSolveResult()
{

}

CLineMatchingSolveResult::~CLineMatchingSolveResult()
{

}

void CLineMatchingSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;

    std::string spectrumName;
    store.GetParam( "spectrumName", spectrumName );

    GetBestResult( store, redshift, merit );

    stream <<  "#Spectrum\tRedshifts\tMatchNum\t"<< std::endl;

    stream  << spectrumName << "\t"
                << redshift << "\t"
                << merit << std::endl;
}

void CLineMatchingSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;

    GetBestResult( store, redshift, merit );
    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << "LineMatchingSolve" << std::endl;
}

Bool CLineMatchingSolveResult::GetBestResult(const CDataStore& store, Float64& redshift, Float64& merit) const
{
    std::string scope = store.GetScope( *this ) + "linematchingsolve.raymatching";
    auto Results =  store.GetGlobalResult(scope.c_str());

    Int32 tmpMatchNum = -1;
    Float64 tmpRedshift = -1.0;

    if(!Results.expired()){
        std::dynamic_pointer_cast<const CRayMatchingResult>( Results.lock() )->GetBestRedshift(tmpRedshift, tmpMatchNum);
    }

    redshift = tmpRedshift;
    merit = tmpMatchNum;
    return true;
}

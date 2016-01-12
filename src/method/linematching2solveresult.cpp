#include <epic/redshift/method/linematching2solveresult.h>
#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <float.h>

using namespace NSEpic;

/**
 * Empty constructor.
 */
CLineMatching2SolveResult::CLineMatching2SolveResult()
{

}

/**
 * Empty destructor.
 */
CLineMatching2SolveResult::~CLineMatching2SolveResult()
{

}

/**
 * Collects from GetBestResult in the store and pretty outputs to stream.
 */
Void CLineMatching2SolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;

    GetBestResult( store, redshift, merit );

    stream << "#Spectrum\tRedshifts\tMatchNum\t"<< std::endl;

    stream << store.GetSpectrumName() << "\t"
	   << redshift << "\t"
	   << merit << std::endl;
}

/**
 * Collects from GetBestResult in the store and pretty outputs to stream, with a "LineMatching2Solve" suffix.
 */
Void CLineMatching2SolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;

    GetBestResult( store, redshift, merit );
    stream << store.GetSpectrumName() << "\t"
	   << redshift << "\t"
	   << merit << "\t"
	   << "LineMatching2Solve" << std::endl;
}

/**
 * Wrapper around CRayMatchingResult::GetBestRedshift.
 */
Bool CLineMatching2SolveResult::GetBestResult( const CDataStore& store, Float64& redshift, Float64& merit ) const
{
    std::string scope = store.GetScope( *this ) + "linematching2solve.raymatching";
    auto Results = store.GetGlobalResult( scope.c_str() );

    Int32 tmpMatchNum = -1;
    Float64 tmpRedshift = -1.0;

    if( !Results.expired() )
      {
        Int32 er = std::dynamic_pointer_cast<const CRayMatchingResult>( Results.lock() )->GetBestRedshift( tmpRedshift, tmpMatchNum );
      }

    redshift = tmpRedshift;
    merit = tmpMatchNum;
    return true;
}


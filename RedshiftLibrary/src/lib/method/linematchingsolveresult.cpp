#include <RedshiftLibrary/method/linematchingsolveresult.h>
#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>

#include <float.h>

using namespace NSEpic;

/**
 * Empty constructor.
 */
CLineMatchingSolveResult::CLineMatchingSolveResult()
{

}

/**
 * Empty destructor.
 */
CLineMatchingSolveResult::~CLineMatchingSolveResult()
{

}

/**
 * Collects from GetBestResult in the store and pretty outputs to stream.
 */
void CLineMatchingSolveResult::Save(std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;

    //TODO review this (commented after removing DataStore from Save and SaveLine)
    /*
      GetBestResult( store, redshift, merit );

    stream << "#Spectrum\tRedshifts\tMatchNum\t"<< std::endl;

    stream << store.GetSpectrumName() << "\t"
	   << redshift << "\t"
	   << merit << std::endl;
    */
}

/**
 * Collects from GetBestResult in the store and pretty outputs to stream, with a "LineMatchingSolve" suffix.
 */
void CLineMatchingSolveResult::SaveLine(std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    //TODO review this (commented after removing DataStore from Save and SaveLine)      
    /*
    GetBestResult( store, redshift, merit );
    stream << store.GetSpectrumName() << "\t"
	   << redshift << "\t"
	   << merit << "\t"
	   << "LineMatchingSolve" << std::endl;
    */
}

/**
 * Wrapper around CRayMatchingResult::GetBestRedshift.
 */
Bool CLineMatchingSolveResult::GetBestResult( const CDataStore& store, Float64& redshift, Float64& merit ) const
{
    std::string scope = store.GetScope( *this ) + "linematchingsolve.raymatching";
    auto Results = store.GetGlobalResult( scope.c_str() );

    Int32 tmpMatchNum = -1;
    Float64 tmpRedshift = -1.0;

    if( !Results.expired() )
      {
        std::dynamic_pointer_cast<const CRayMatchingResult>( Results.lock() )->GetBestRedshift( tmpRedshift, tmpMatchNum );
      }

    redshift = tmpRedshift;
    merit = tmpMatchNum;
    return true;
}


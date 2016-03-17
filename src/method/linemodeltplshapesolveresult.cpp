#include <epic/redshift/method/linemodeltplshapesolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/linemodelresult.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CLineModelTplshapeSolveResult::CLineModelTplshapeSolveResult()
{

}

/**
 * \brief Empty destructor.
 **/
CLineModelTplshapeSolveResult::~CLineModelTplshapeSolveResult()
{

}

/**
 * \brief Outputs to the output stream the values for redshift and merit and template name of the best redshift obtained.
 **/
Void CLineModelTplshapeSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit );

    stream <<  "#Redshifts\tMerit\tTemplate"<< std::endl;
    stream << redshift << "\t"
	   << merit << "\t"
	   << tplName << std::endl;
}

/**
 * \brief Prints into the output stream the redshift, merit and template name for the best redshift obtained.
 **/
Void CLineModelTplshapeSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit );

    stream  << store.GetSpectrumName() << "\t"
	    << redshift << "\t"
	    << merit << "\t"
	    << tplName << "\t"
        << "LineModeltplshapeSolve" << std::endl;
}

/**
 * \brief Searches all the results for the first one with the lowest value of merit, and uses it to update the argument references.
 * Construct the scope string for the Linemodel results.
 * Set temporary variables with maximum merit and 0 redshift.
 * If the results are not expired, lock the mutex, and for each entry in the result chisquare:
 *   if this entry is less than the temporary merit, update the temporary merit and redshift values with the result stored in this entry.
 * Set the redshift and merit referenced arguments with the best values found.
 **/
Bool CLineModelTplshapeSolveResult::GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit ) const
{
    std::string scope = store.GetScope( *this ) + "linemodeltplshapesolve.linemodel";
    auto results = store.GetGlobalResult( scope.c_str() );

    Float64 tmpMerit = DBL_MAX;
    Float64 tmpRedshift = 0.0;

    if( !results.expired() )
      {
        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );
        for( Int32 i=0; i<lineModelResult->ChiSquare.size(); i++ )
	  {
            if( lineModelResult->ChiSquare[i] < tmpMerit )
	      {
                tmpMerit = lineModelResult->ChiSquare[i];
                tmpRedshift = lineModelResult->Redshifts[i];
	      }
	  }
      }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
}

/**
 * \brief Searches all the results for the first one with the largest value for the LogArea, and uses it to update the argument references.
 * Construct the scope string for the Linemodel results.
 * Set temporary variables with minimum merit and 0 redshift.
 * If the results are not expired, lock the mutex, and for each entry in the result LogArea:
 *   if this entry is larger than the temporary merit, update the temporary merit and redshift values with the results stored in this entry.
 * Set the redshift and merit referenced arguments with the best values found.
 **/
Bool CLineModelTplshapeSolveResult::GetBestRedshiftLogArea( const CDataStore& store, Float64& redshift, Float64& merit ) const
{
    std::string scope = store.GetScope( *this ) + "linemodelsolve.linemodel";
    auto results = store.GetGlobalResult( scope.c_str() );

    Float64 tmpMerit = -DBL_MAX ;
    Float64 tmpRedshift = 0.0;

    if(!results.expired())
      {
        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );
        for( Int32 i=0; i<lineModelResult->LogArea.size(); i++ )
	  {
            if( lineModelResult->LogArea[i] > tmpMerit )
            {
                tmpMerit = lineModelResult->LogArea[i];
                tmpRedshift = lineModelResult->LogAreaCorrectedExtrema[i];
            }
	  }
      }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
}




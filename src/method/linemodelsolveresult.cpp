#include <epic/redshift/method/linemodelsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/linemodelresult.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CLineModelSolveResult::CLineModelSolveResult()
{

}

/**
 * \brief Empty destructor.
 **/
CLineModelSolveResult::~CLineModelSolveResult()
{

}

/**
 * \brief Outputs to the output stream the values for redshift and merit and template name of the best redshift obtained.
 **/
Void CLineModelSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    Float64 sigma;

    //GetBestRedshiftWithStrongELSnrPrior( store, redshift, merit );
    GetBestRedshift( store, redshift, merit, sigma );

    stream <<  "#Redshifts\tMerit\tTemplate"<< std::endl;
    stream << redshift << "\t"
	   << merit << "\t"
       << tplName << "\t"
       << "LineModelSolve" << "\t"
       << sigma << std::endl;
}

/**
 * \brief Prints into the output stream the redshift, merit and template name for the best redshift obtained.
 **/
Void CLineModelSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    Float64 sigma;

    //GetBestRedshiftWithStrongELSnrPrior( store, redshift, merit );
    GetBestRedshift( store, redshift, merit, sigma );

    stream  << store.GetSpectrumName() << "\t"
	    << redshift << "\t"
	    << merit << "\t"
	    << tplName << "\t"
        << "LineModelSolve" << "\t"
        << sigma << std::endl;
}

/**
 * \brief Searches all the results for the first one with the lowest value of merit, and uses it to update the argument references.
 * Construct the scope string for the Linemodel results.
 * Set temporary variables with maximum merit and 0 redshift.
 * If the results are not expired, lock the mutex, and for each entry in the result chisquare:
 *   if this entry is less than the temporary merit, update the temporary merit and redshift values with the result stored in this entry.
 * Set the redshift and merit referenced arguments with the best values found.
 **/
Bool CLineModelSolveResult::GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, Float64& sigma ) const
{
    std::string scope = store.GetScope( *this ) + "linemodelsolve.linemodel";
    auto results = store.GetGlobalResult( scope.c_str() );

    Float64 tmpMerit = DBL_MAX;
    Float64 tmpRedshift = 0.0;
    Float64 tmpSigma = -1.0;

    if( !results.expired() )
    {
        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );
        for( Int32 i=0; i<lineModelResult->Extrema.size(); i++ )
        {
            Float64 merit = lineModelResult->GetExtremaMerit(i);
            if( merit < tmpMerit )
            {
                tmpMerit = merit;
                tmpRedshift = lineModelResult->Extrema[i];
                tmpSigma = lineModelResult->DeltaZ[i];
            }
        }
    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    sigma = tmpSigma;
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
Bool CLineModelSolveResult::GetBestRedshiftLogArea( const CDataStore& store, Float64& redshift, Float64& merit ) const
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

/**
 * \brief Searches all the extrema results for the lowest, using the StrongELSnrPrior to correct the initial Chi2 value
 **/
Bool CLineModelSolveResult::GetBestRedshiftWithStrongELSnrPrior( const CDataStore& store, Float64& redshift, Float64& merit ) const
{
    std::string scope = store.GetScope( *this ) + "linemodelsolve.linemodel";
    auto results = store.GetGlobalResult( scope.c_str() );

    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;

    if(!results.expired())
    {
        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );
        for( Int32 i=0; i<lineModelResult->Extrema.size(); i++ )
        {
            Float64 coeff =10.0;
            Float64 correctedMerit = lineModelResult->GetExtremaMerit(i)-coeff*lineModelResult->StrongELSNR[i];
            if( correctedMerit < tmpMerit )
            {
                tmpMerit = correctedMerit;
                tmpRedshift = lineModelResult->Extrema[i];
            }
        }
    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
}




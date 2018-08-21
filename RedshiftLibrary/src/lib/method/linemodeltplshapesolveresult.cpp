#include <RedshiftLibrary/method/linemodeltplshapesolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/linemodelresult.h>
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
void CLineModelTplshapeSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit, tplName);

    stream <<  "#Redshifts\tMerit\tTemplateShape"<< std::endl;
    stream << redshift << "\t"
	   << merit << "\t"
	   << tplName << std::endl;



    stream << std::endl;
    stream << std::endl;
    std::string detailStr;
    GetBestRedshiftPerTemplateString( store, detailStr);

    stream << detailStr.c_str();
}

/**
 * \brief Prints into the output stream the redshift, merit and template name for the best redshift obtained.
 **/
void CLineModelTplshapeSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshift( store, redshift, merit, tplName);

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
Bool CLineModelTplshapeSolveResult::GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit , std::string& tplName) const
{
    std::string scope = store.GetScope( *this ) + "linemodeltplshapesolve.linemodel";
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CLineModelResult>( (*it).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->ChiSquare[i] < tmpMerit && meritResult->Status[i] == COperator::nStatus_OK )
            {
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpTplName = (*it).first;
            }
        }
    }

    if( tmpMerit < DBL_MAX )
    {
        redshift = tmpRedshift;
        merit = tmpMerit;
        tplName = tmpTplName;
        return true;
    }

    return false;
}

Bool CLineModelTplshapeSolveResult::GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const
{
    std::string scope = store.GetScope( *this ) + "linemodeltplshapesolve.linemodel";
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());


    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        Float64 tmpMerit = DBL_MAX ;
        Float64 tmpRedshift = 0.0;
        std::string tmpTplName;

        auto meritResult = std::dynamic_pointer_cast<const CLineModelResult>( (*it).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->ChiSquare[i] < tmpMerit && meritResult->Status[i] == COperator::nStatus_OK )
            {
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpTplName = (*it).first;
            }
        }



        if( tmpMerit < DBL_MAX )
        {
            char tmpChar[256];
            sprintf(tmpChar, "%f\t%f\t%s\n", tmpRedshift, tmpMerit, tmpTplName.c_str());
            output.append(tmpChar);
        }else{
            char tmpChar[256];
            sprintf(tmpChar, "-1\t-1\t%s\n", tmpTplName.c_str());
            output.append(tmpChar);
        }
    }


    return true;

}



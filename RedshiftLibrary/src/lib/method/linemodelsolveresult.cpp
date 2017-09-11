#include <RedshiftLibrary/method/linemodelsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/linemodelresult.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>


using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CLineModelSolveResult::CLineModelSolveResult()
{
    m_ReliabilityLabel = "-1";
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
    std::string tplName="-1";
    Float64 sigma;

    //GetBestRedshiftWithStrongELSnrPrior( store, redshift, merit );
    if(m_bestRedshiftMethod==0)
    {
        GetBestRedshift( store, redshift, merit, sigma );
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from chi2 extrema: z=%f", redshift);
    }else if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit, sigma );
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from PDF: z=%f", redshift);
    }else{
        Log.LogError( "Linemodelsolve-result: can't parse best redshift estimation method");
    }


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
    std::string tplName="-1";
    Float64 sigma;


    //GetBestRedshiftWithStrongELSnrPrior( store, redshift, merit );
    if(m_bestRedshiftMethod==0)
    {
        GetBestRedshift( store, redshift, merit, sigma );
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from chi2 extrema: z=%f", redshift);;
    }else if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit, sigma );
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from PDF: z=%f", redshift);
    }else{
        Log.LogError( "Linemodelsolve-result: can't parse best redshift estimation method");
    }

    stream  << store.GetSpectrumName() << "\t"
        << store.GetProcessingID() << "\t"
	    << redshift << "\t"
	    << merit << "\t"
	    << tplName << "\t"
        << "LineModelSolve" << "\t"
        << sigma << "\t"
        << m_ReliabilityLabel << std::endl;
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
        for( Int32 i=0; i<lineModelResult->ExtremaResult.Extrema.size(); i++ )
        {
            Float64 merit = lineModelResult->GetExtremaMerit(i);
            if( merit < tmpMerit )
            {
                tmpMerit = merit;
                //tmpRedshift = lineModelResult->Extrema[i];
                tmpRedshift = lineModelResult->ExtremaResult.ExtremaLastPass[i];
                //tmpRedshift = lineModelResult->lmfitPass[i];
                tmpSigma = lineModelResult->ExtremaResult.DeltaZ[i];
            }
        }
    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    sigma = tmpSigma;
    return true;
}

/**
 * @brief CLineModelSolveResult::isPdfValid
 * @return
 */
Bool CLineModelSolveResult::isPdfValid(const CDataStore& store) const
{
    std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
    auto results_pdf =  store.GetGlobalResult( scope_res.c_str() );
    auto logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results_pdf.lock() );

    if(!logzpdf1d)
    {
        return false;
    }

    if(logzpdf1d->Redshifts.size()<2)
    {
        return false;
    }

    return true;
}


/**
 * \brief Searches the best_z = argmax(pdf)
 * output: redshift = argmax(pdf)
 * output: merit = chi2(redshift)
 * output: sigma = deltaz(redshift)
 *
 **/
Bool CLineModelSolveResult::GetBestRedshiftFromPdf( const CDataStore& store, Float64& redshift, Float64& merit, Float64& sigma ) const
{
    std::string scope = store.GetScope( *this ) + "linemodelsolve.linemodel";
    auto results_chi2 = store.GetGlobalResult( scope.c_str() );

    std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
    auto results_pdf =  store.GetGlobalResult( scope_res.c_str() );
    auto logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results_pdf.lock() );

    if(!logzpdf1d)
    {
        Log.LogError( "GetBestRedshiftFromPdf: no pdf results retrieved from scope: %s", scope_res.c_str());
        return false;
    }



    Float64 tmpProbaLog = -DBL_MAX;
    Float64 tmpMerit = DBL_MAX;
    Float64 tmpRedshift = 0.0;
    Float64 tmpSigma = -1.0;

    if( !results_chi2.expired() )
    {
        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results_chi2.lock() );


        if(logzpdf1d->Redshifts.size() != lineModelResult->Redshifts.size())
        {
            Log.LogError( "GetBestRedshiftFromPdf: pdf samplecount != chisquare samplecount");
            return false;
        }

        for( Int32 i=0; i<lineModelResult->ExtremaResult.Extrema.size(); i++ )
        {
            UInt32 solIdx = lineModelResult->GetExtremaIndex(i);
            if(solIdx<0 || solIdx>=logzpdf1d->valProbaLog.size())
            {
                Log.LogError( "GetBestRedshiftFromPdf: pdf proba value not found for extremumIndex = %d", i);
                return false;
            }

            Float64 probaLog = logzpdf1d->valProbaLog[solIdx];

            Float64 merit = lineModelResult->GetExtremaMerit(i);
            //if( merit < tmpMerit )
            if(probaLog>tmpProbaLog)
            {
                tmpProbaLog = probaLog;
                tmpMerit = merit;
                //tmpRedshift = lineModelResult->Extrema[i];
                tmpRedshift = lineModelResult->ExtremaResult.ExtremaLastPass[i];
                tmpSigma = lineModelResult->ExtremaResult.DeltaZ[i];
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
        for( Int32 i=0; i<lineModelResult->ExtremaResult.LogArea.size(); i++ )
      {
            if( lineModelResult->ExtremaResult.LogArea[i] > tmpMerit )
            {
                tmpMerit = lineModelResult->ExtremaResult.LogArea[i];
                tmpRedshift = lineModelResult->ExtremaResult.LogAreaCorrectedExtrema[i];
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
        for( Int32 i=0; i<lineModelResult->ExtremaResult.Extrema.size(); i++ )
        {
            Float64 coeff =10.0;
            Float64 correctedMerit = lineModelResult->GetExtremaMerit(i)-coeff*lineModelResult->ExtremaResult.StrongELSNR[i];
            if( correctedMerit < tmpMerit )
            {
                tmpMerit = correctedMerit;
                tmpRedshift = lineModelResult->ExtremaResult.Extrema[i];
            }
        }
    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
}

Void CLineModelSolveResult::SetReliabilityLabel( std::string lbl )
{
    m_ReliabilityLabel = lbl;
}




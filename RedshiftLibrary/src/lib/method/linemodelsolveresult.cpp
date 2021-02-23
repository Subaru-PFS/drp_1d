#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/linemodelresult.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/operator/pdfz.h>

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CLineModelSolveResult::CLineModelSolveResult(const std::shared_ptr<const CLineModelExtremaResult> & ExtremaResult,
                                             const std::string & opt_pdfcombination,
                                             Float64 evidence):
    CSolveResult( ExtremaResult, opt_pdfcombination, evidence),
    ExtremaResult(ExtremaResult),
    tplratioName(ExtremaResult->FittedTplratioName[0]),
    tplcontinuumName(ExtremaResult->FittedTplName[0]),
    sigma(ExtremaResult->DeltaZ(0)),
    snrHa(ExtremaResult->snrHa[0]),
    lfHa(ExtremaResult->lfHa[0]),
    snrOII(ExtremaResult->snrOII[0]),
    lfOII(ExtremaResult->lfOII[0])
{}

/**
 * \brief Empty destructor.
 **/
CLineModelSolveResult::~CLineModelSolveResult()
{


}
/*
void CLineModelSolveResult::preSave(const CDataStore& store)
{
   //GetBestRedshiftWithStrongELSnrPrior( store, redshift, merit );
    if(m_bestRedshiftMethod==0)
    {
        GetBestRedshift( m_redshift, m_merit, sigma, snrHa, lfHa, snrOII, lfOII );
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from chi2 extrema: z=%f", m_redshift);
    }
    else if(m_bestRedshiftMethod==2)
    {
        if(GetBestRedshiftFromPdf( m_redshift, m_merit, sigma, snrHa, lfHa, snrOII, lfOII, tplratioName, tplcontinuumName ))
        {
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from PDF: z=%f", m_redshift);
        Log.LogInfo( "Linemodelsolve-result: extracted best model-tplratio=%s", tplratioName.c_str());
        Log.LogInfo( "Linemodelsolve-result: extracted best model-tplcontinuum=%s", tplcontinuumName.c_str());
        }
        else Log.LogError( "Linemodelsolve-result: can't get best redshift From Pdf");
    }
    else{
        Log.LogError( "Linemodelsolve-result: can't parse best redshift estimation method");
    }

}*/

/**
 * \brief Outputs to the output stream the values for redshift and merit and template name of the best redshift obtained.
 **/
void CLineModelSolveResult::Save(std::ostream& stream ) const
{
    stream <<  "#Redshifts\tMerit\tTemplateRatio\tTemplateContinuum\tmethod\tsigma"<< std::endl;
    stream << m_redshift << "\t"
       << m_merit << "\t"
       << tplratioName << "\t"
       << tplcontinuumName << "\t"
       << "LineModelSolve" << "\t"
       << sigma << std::endl;
}

/**
 * \brief Prints into the output stream the redshift, merit and template name for the best redshift obtained.
 **/
void CLineModelSolveResult::SaveLine( std::ostream& stream ) const
{/*
    Float64 redshift;
    Float64 merit;
    std::string tplratioName="-1";
    std::string tplcontinuumName="-1";
    Float64 sigma;
    Float64 snrHa;
    Float64 lfHa=-1;
    Float64 snrOII;
    Float64 lfOII=-1;


    //GetBestRedshiftWithStrongELSnrPrior( store, redshift, merit );
    if(m_bestRedshiftMethod==0)
    {
    //TODO review this (commented after removing DataStore from Save and SaveLine)      
        //        GetBestRedshift( store, redshift, merit, sigma, snrHa, lfHa, snrOII, lfOII );
        //Log.LogInfo( "Linemodelsolve-result: extracting best redshift from chi2 extrema: z=%f", redshift);;
    }else if(m_bestRedshiftMethod==2)
    {
          //TODO review this (commented after removing DataStore from Save and SaveLine)
      //        GetBestRedshiftFromPdf( store, redshift, merit, sigma, snrHa, lfHa, snrOII, lfOII, tplratioName, tplcontinuumName);
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from PDF: z=%f", redshift);
    }else{
        //Log.LogError( "Linemodelsolve-result: can't parse best redshift estimation method");
    }
    //TODO review this (commented after removing DataStore from Save and SaveLine)      
    /*
    stream  << store.GetSpectrumName() << "\t"
        << store.GetProcessingID() << "\t"
	    << m_redshift << "\t"
	    << m_merit << "\t"
        << tplratioName << "\t"
        << "LineModelSolve" << "\t"
        << sigma << "\t"
        << m_ReliabilityLabel << "\t"
        << snrHa << "\t"
        << lfHa << "\t"
        << snrOII << "\t"
        << lfOII << "\t"
        << m_TypeLabel << std::endl;
    */
}


/**
 * \brief Searches all the results for the first one with the lowest value of merit, and uses it to update the argument references.
 * Construct the scope string for the Linemodel results.
 * Set temporary variables with maximum merit and 0 redshift.
 * If the results are not expired, lock the mutex, and for each entry in the result chisquare:
 *   if this entry is less than the temporary merit, update the temporary merit and redshift values with the result stored in this entry.
 * Set the redshift and merit referenced arguments with the best values found.
 **/
/*
Bool CLineModelSolveResult::GetBestRedshift(Float64& redshift,
                                            Float64& merit,
                                            Float64& sigma,
                                            Float64& snrHa,
                                            Float64& lfHa ,
                                            Float64 &snrOII,
                                            Float64 &lfOII) const
{
    Float64 tmpMerit = DBL_MAX;
    Float64 tmpRedshift = 0.0;
    Float64 tmpSigma = -1.0;
    Float64 tmpSnrHa = -1.0;
    Float64 tmpLFHa = -1.0;
    Float64 tmpSnrOII = -1.0;
    Float64 tmpLFOII = -1.0;
    
    for( Int32 i=0; i<ExtremaResult->size(); i++ )
    {
        Float64 merit = ExtremaResult->ValProba(i); 
        if( merit > tmpMerit )
        {
            tmpMerit = merit;
            tmpRedshift = ExtremaResult->Redshift(i);
            //tmpRedshift = ExtremaResult->Redshift_lmfit[i];
            tmpSigma = ExtremaResult->DeltaZ(i);
            tmpSnrHa = ExtremaResult->snrHa[i];
            tmpLFHa = ExtremaResult->lfHa[i];
            tmpSnrOII = ExtremaResult->snrOII[i];
            tmpLFOII = ExtremaResult->lfOII[i];
        }
    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    sigma = tmpSigma;
    snrHa = tmpSnrHa;
    lfHa = tmpLFHa;
    snrOII = tmpSnrOII;
    lfOII = tmpLFOII;
    return true;
}*/

/**
 * Simply reading from datastore info related to the best Candidate
*/
/*
Bool CLineModelSolveResult::GetBestRedshiftFromPdf( Float64& redshift,
                                                    Float64& probaLog,
                                                    Float64& sigma,
                                                    Float64& snrHa,
                                                    Float64& lfHa ,
                                                    Float64 &snrOII,
                                                    Float64 &lfOII,
                                                    std::string& modelTplratio,
                                                    std::string& modelTplContinuum ) const
{    
    Int32 bestIdx = 0;

    const TCandidateZ & bestZCandidate = ExtremaResult->Candidates[bestIdx].second;

    redshift = bestZCandidate.Redshift;
    probaLog = bestZCandidate.ValSumProba;
    sigma =   bestZCandidate.Deltaz;
    snrHa = ExtremaResult->snrHa[bestIdx];
    lfHa = ExtremaResult->lfHa[bestIdx];
    snrOII = ExtremaResult->snrOII[bestIdx];
    lfOII = ExtremaResult->lfOII[bestIdx];
    modelTplratio = ExtremaResult->FittedTplratioName[bestIdx];
    modelTplContinuum = ExtremaResult->FittedTplName[bestIdx];
    return true;
}*/


/*Int32 CLineModelSolveResult::GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) 
{
    //check if the pdf is stellar/galaxy or else
    std::string outputPdfRelDir  = "zPDF"; //default value set to galaxy zPDF folder
    std::string scope = store.GetScope( *this );
    std::size_t foundstr_stellar = scope.find("stellarsolve");
    std::size_t foundstr_qso = scope.find("qsosolve");
    if (foundstr_stellar!=std::string::npos){
        outputPdfRelDir = "stellar_zPDF";
    }else if (foundstr_qso!=std::string::npos){
        outputPdfRelDir = "qso_zPDF";
    }


    std::string scope_res = outputPdfRelDir+"/logposterior.logMargP_Z_data";
    auto results_pdf =  store.GetGlobalResult( scope_res.c_str() );
    std::shared_ptr<const CPdfMargZLogResult> logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results_pdf.lock() );

    if(!logzpdf1d)
    {
        Log.LogError( "GetEvidenceFromPdf: no pdf results retrieved from scope: %s", scope_res.c_str());
        return -1;
    }

    evidence = logzpdf1d->valEvidenceLog;
    return 0;
}*/

/**
 * \brief Searches all the results for the first one with the largest value for the LogArea, and uses it to update the argument references.
 * Construct the scope string for the Linemodel results.
 * Set temporary variables with minimum merit and 0 redshift.
 * If the results are not expired, lock the mutex, and for each entry in the result LogArea:
 *   if this entry is larger than the temporary merit, update the temporary merit and redshift values with the results stored in this entry.
 * Set the redshift and merit referenced arguments with the best values found.
 **/
/*
Bool CLineModelSolveResult::GetBestRedshiftLogArea( Float64& redshift, Float64& merit ) const
{
    Float64 tmpMerit = -DBL_MAX ;
    Float64 tmpRedshift = 0.0;

    for( Int32 i=0; i<ExtremaResult->LogArea.size(); i++ )
      {
            if( ExtremaResult->LogArea[i] > tmpMerit )
            {
                tmpMerit = ExtremaResult->LogArea[i];
                tmpRedshift = ExtremaResult->LogAreaCorrectedExtrema[i];
            }
      }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
}*/

void CLineModelSolveResult::getData(const std::string& name, Float64& v) const
{
  if (name.compare("snrHa") == 0)  v = snrHa;
  else if (name.compare("lfHa") == 0)  v = lfHa;
  else if (name.compare("snrOII") == 0)  v = snrOII;
  else if (name.compare("lfOII") == 0)  v = lfOII;
}

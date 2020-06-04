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
#include <RedshiftLibrary/statistics/pdfz.h>

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
void CLineModelSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplratioName="-1";
    std::string tplcontinuumName="-1";
    Float64 sigma;
    Float64 snrHa=-1.0;
    Float64 lfHa=-1.0;
    Float64 snrOII=-1.0;
    Float64 lfOII=-1.0;

    //GetBestRedshiftWithStrongELSnrPrior( store, redshift, merit );
    if(m_bestRedshiftMethod==0)
    {
        GetBestRedshift( store, redshift, merit, sigma, snrHa, lfHa, snrOII, lfOII );
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from chi2 extrema: z=%f", redshift);
    }else if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit, sigma, snrHa, lfHa, snrOII, lfOII, tplratioName, tplcontinuumName );
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from PDF: z=%f", redshift);
        Log.LogInfo( "Linemodelsolve-result: extracted best model-tplratio=%s", tplratioName.c_str());
        Log.LogInfo( "Linemodelsolve-result: extracted best model-tplcontinuum=%s", tplcontinuumName.c_str());
    }else{
        Log.LogError( "Linemodelsolve-result: can't parse best redshift estimation method");
    }


    stream <<  "#Redshifts\tMerit\tTemplateRatio\tTemplateContinuum\tmethod\tsigma"<< std::endl;
    stream << redshift << "\t"
       << merit << "\t"
       << tplratioName << "\t"
       << tplcontinuumName << "\t"
       << "LineModelSolve" << "\t"
       << sigma << std::endl;
}

/**
 * \brief Prints into the output stream the redshift, merit and template name for the best redshift obtained.
 **/
void CLineModelSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
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
        GetBestRedshift( store, redshift, merit, sigma, snrHa, lfHa, snrOII, lfOII );
        //Log.LogInfo( "Linemodelsolve-result: extracting best redshift from chi2 extrema: z=%f", redshift);;
    }else if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit, sigma, snrHa, lfHa, snrOII, lfOII, tplratioName, tplcontinuumName);
        Log.LogInfo( "Linemodelsolve-result: extracting best redshift from PDF: z=%f", redshift);
    }else{
        //Log.LogError( "Linemodelsolve-result: can't parse best redshift estimation method");
    }


    stream  << store.GetSpectrumName() << "\t"
        << store.GetProcessingID() << "\t"
	    << redshift << "\t"
	    << merit << "\t"
        << tplratioName << "\t"
        << "LineModelSolve" << "\t"
        << sigma << "\t"
        << m_ReliabilityLabel << "\t"
        << snrHa << "\t"
        << lfHa << "\t"
        << snrOII << "\t"
        << lfOII << "\t"
        << m_TypeLabel << std::endl;
}

/**
 * \brief Searches all the results for the first one with the lowest value of merit, and uses it to update the argument references.
 * Construct the scope string for the Linemodel results.
 * Set temporary variables with maximum merit and 0 redshift.
 * If the results are not expired, lock the mutex, and for each entry in the result chisquare:
 *   if this entry is less than the temporary merit, update the temporary merit and redshift values with the result stored in this entry.
 * Set the redshift and merit referenced arguments with the best values found.
 **/
Bool CLineModelSolveResult::GetBestRedshift(const CDataStore& store,
                                            Float64& redshift,
                                            Float64& merit,
                                            Float64& sigma,
                                            Float64& snrHa,
                                            Float64& lfHa ,
                                            Float64 &snrOII,
                                            Float64 &lfOII) const
{
    std::string scope = store.GetScope( *this ) + "linemodelsolve.linemodel";
    auto results = store.GetGlobalResult( scope.c_str() );

    Float64 tmpMerit = DBL_MAX;
    Float64 tmpRedshift = 0.0;
    Float64 tmpSigma = -1.0;
    Float64 tmpSnrHa = -1.0;
    Float64 tmpLFHa = -1.0;
    Float64 tmpSnrOII = -1.0;
    Float64 tmpLFOII = -1.0;

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
                tmpSnrHa = lineModelResult->ExtremaResult.snrHa[i];
                tmpLFHa = lineModelResult->ExtremaResult.lfHa[i];
                tmpSnrOII = lineModelResult->ExtremaResult.snrOII[i];
                tmpLFOII = lineModelResult->ExtremaResult.lfOII[i];
            }
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
}

Bool CLineModelSolveResult::GetBestRedshiftsFromPdf(const CDataStore& store,
                                                    TFloat64List& redshift) const
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

    if( !results_chi2.expired() )
    {
        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results_chi2.lock() );

        if(logzpdf1d->Redshifts.size() != lineModelResult->Redshifts.size())
        {
            Log.LogError( "GetBestRedshiftFromPdf: pdf samplecount != chisquare samplecount");
            return false;
        }
        Float64 Fullwidth = 6e-3;//should be replaced with deltaz?
        Int32 method = 1; //0=maxpdf, 1=direct integration on peaks only
        for( Int32 i=0; i<lineModelResult->ExtremaResult.Extrema.size(); i++ )
        {
            Float64 tmpIntgProba = -DBL_MAX;
            Float64 tmpRedshift = 0.0; 
            /*//TODO: below commented code should replace executing code once the new finder is ready; find() arguments shd also be updated
            TPointList extremumList;
            Int32 s = lineModelResult->ExtremaResult.ExtremaExtendedRedshifts[i].size();
            TFloat64Range redshiftsRange = TFloat64Range( lineModelResult->ExtremaResult.ExtremaExtendedRedshifts[i][0], 
                                                           lineModelResult->ExtremaResult.ExtremaExtendedRedshifts[i][s-1]);
            //call Find on each secondpass range and retrieve the best 10 peaks?
            CExtremum extremum(redshiftsRange, 10, false, 1);
            extremum.Find(logzpdf1d->Redshifts, logzpdf1d->valProbaLog, extremumList);
            if(method == 0){
                redshift.push_back(extremumList[0].X);
                continue;
            }
            for(Int32 j = 0; j < extremumList.size(); j++){
                CPdfz pdfz;
                Float64 flux_integral = -1;
                flux_integral = pdfz.getCandidateSumTrapez( logzpdf1d->Redshifts, logzpdf1d->valProbaLog, extremumList[j].X, Fullwidth);
                if(flux_integral>tmpIntgProba){
                        tmpRedshift = extremumList[j].X;
                        tmpIntgProba = extremumList[j].Y;
                }
            }*/
            Float64 tmpProbaLog = -DBL_MAX;
            for(Int32 kval=0; kval<lineModelResult->ExtremaResult.ExtremaExtendedRedshifts[i].size(); kval++)
            {
                Float64 zInCandidateRange = lineModelResult->ExtremaResult.ExtremaExtendedRedshifts[i][kval];
                UInt32 solIdx = logzpdf1d->getIndex(zInCandidateRange);
                if(solIdx<0 || solIdx>=logzpdf1d->valProbaLog.size())
                {
                    Log.LogError( "GetBestRedshiftFromPdf: pdf proba value not found for extremumIndex = %d", i);
                    return false;
                }

                Float64 probaLog = logzpdf1d->valProbaLog[solIdx];
                Log.LogDebug( "GetBestRedshiftFromPdf: z=%f : probalog = %f", zInCandidateRange, probaLog);
                
                if(method == 0){
                    if(probaLog>tmpProbaLog){
                        tmpRedshift = zInCandidateRange;
                        tmpProbaLog = probaLog;
                    }
                }    
                if(method == 1){//max integrated proba but only on peaks in this range
                    CPdfz pdfz;
                    Float64 flux_integral = -1;
                    Float64 prev, next;
                    prev = logzpdf1d->valProbaLog[solIdx-1];
                    next = logzpdf1d->valProbaLog[solIdx+1];
                    if((probaLog > prev&& probaLog > next) ||
                        (solIdx == 0 && probaLog>next) ||
                        (solIdx == logzpdf1d->valProbaLog.size()-1 && probaLog > prev)){
                        //if current value is a peak
                        flux_integral = pdfz.getCandidateSumTrapez( logzpdf1d->Redshifts, logzpdf1d->valProbaLog, zInCandidateRange, Fullwidth);
                        if(flux_integral>tmpIntgProba){
                            tmpRedshift = zInCandidateRange;
                            tmpIntgProba = probaLog;
                        }
                    }
                    else 
                    {
                        continue; //it doesnt work to compute here the pdfz
                    }
                } 
            }
            redshift.push_back(tmpRedshift);
        }

    }

    return true;
}
/**
 * Simply reading from datastore info related to the best Candidate
*/
Bool CLineModelSolveResult::GetBestRedshiftFromPdf(const CDataStore& store,
                                            Float64& redshift,
                                            Float64& probaLog,
                                            Float64& sigma,
                                            Float64& snrHa,
                                            Float64& lfHa ,
                                            Float64 &snrOII,
                                            Float64 &lfOII,
                                            std::string& modelTplratio,
                                            std::string& modelTplContinuum ) const
{
    std::string scope = store.GetScope( *this ) + "linemodelsolve.linemodel";
    auto results = store.GetGlobalResult( scope.c_str() );
    auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );
    //ideally ExtremaPDF should be saved in resultstore as part of lineModelResult object! 
    
    scope = store.GetScope( *this ) + "candidatesresult";
    auto res = store.GetGlobalResult( scope.c_str() );
    auto candResults = std::dynamic_pointer_cast<const CPdfCandidateszResult>( res.lock());
    TFloat64List ExtremaPDF = candResults->ValSumProba;
    //reading from candResults cause deltaz is computed correctly there and what is saved in datastore is bad
    TFloat64List ExtremaDeltaz = candResults->Deltaz;
    TFloat64List Extrema = candResults->Redshifts;
    Int32 bestIdx = candResults->Rank[0];

    if(bestIdx>=Extrema.size() || !Extrema.size()){
       Log.LogError("CLineModelSolveResult::GetBestRedshiftFromPdf: Can't access best redshift ");
       throw runtime_error("Can't access best redshift");
    }
    if(results.expired())
        return false;
    //is not possible, we are reading values from datastore that we update in pdfzcandidatesresult!!!
    redshift = Extrema[bestIdx];
    probaLog = ExtremaPDF[bestIdx];
    sigma = ExtremaDeltaz[bestIdx];
    //not sure that below values are correct!
    snrHa = lineModelResult->ExtremaResult.snrHa[bestIdx];
    lfHa = lineModelResult->ExtremaResult.lfHa[bestIdx];
    snrOII = lineModelResult->ExtremaResult.snrOII[bestIdx];
    lfOII = lineModelResult->ExtremaResult.lfOII[bestIdx];
    modelTplratio = lineModelResult->ExtremaResult.FittedTplshapeName[bestIdx];
    return true;
}


Int32 CLineModelSolveResult::GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
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

Bool CLineModelSolveResult::GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates) const
{
    Log.LogDebug( "CLineModelSolveResult:GetRedshiftCandidates" );
    redshiftcandidates.clear();
    std::string scope = store.GetScope( *this ) + "linemodelsolve.linemodel";
    auto results = store.GetGlobalResult( scope.c_str() );

    if(!results.expired())
    {
        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );
        for( Int32 i=0; i<lineModelResult->ExtremaResult.Extrema.size(); i++ )
        {
            Float64 tmpRedshift = lineModelResult->ExtremaResult.Extrema[i];
            redshiftcandidates.push_back(tmpRedshift);
        }
    }else{
        return false;
    }

    return true;
}


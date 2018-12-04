#include <RedshiftLibrary/method/chisquare2solveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>
#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/extremum/extremum.h>

using namespace NSEpic;

CChisquare2SolveResult::CChisquare2SolveResult()
{
    m_type = nType_raw;
}

CChisquare2SolveResult::~CChisquare2SolveResult()
{

}

void CChisquare2SolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    Float64 evidence;
    std::string tplName;
    Float64 amplitude;
    Float64 dustCoeff;
    Int32 meiksinIdx;

    GetBestRedshift( store, redshift, merit, tplName, amplitude, dustCoeff, meiksinIdx );

    stream <<  "#Redshifts\tMerit\tTemplate\tampl\tdustcoeff\tmeiksinidx"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << amplitude << "\t"
                << std::setprecision(4) << dustCoeff << "\t"
                << meiksinIdx << std::endl;

    stream <<  "#Redshifts\tprobaLog\tevidenceLog\tModel"<< std::endl;
    if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit, evidence );
        Log.LogInfo( "Chisquare2solve-result: extracting best redshift from PDF: z=%f", redshift);
        GetBestModel(store, redshift, tplName);
        Log.LogInfo( "Chisquare2solve-result: extracted best model: model=%s", tplName.c_str());

        stream  << redshift << "\t"
                << merit << "\t"
                << evidence << "\t"
                << tplName << std::endl;
    }else{
        stream <<  "-1\t-1\t-1"<< std::endl;
    }


    stream << std::endl;
    stream << std::endl;
    std::string detailStr;
    GetBestRedshiftPerTemplateString( store, detailStr);

    stream << detailStr.c_str();
}

Bool CChisquare2SolveResult::GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const
{
    std::string scopeStr;
    if(m_type == nType_raw){
        scopeStr = "chisquare";
    }else if(m_type == nType_all){
        scopeStr = "chisquare";
    }else if(m_type == nType_noContinuum){
        scopeStr = "chisquare_nocontinuum";
    }else if(m_type == nType_continuumOnly){
        scopeStr = "chisquare_continuum";
    }
    std::string scope = store.GetScope( *this ) + "chisquare2solve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());


    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        Float64 tmpMerit = DBL_MAX ;
        Float64 tmpRedshift = 0.0;
        std::string tmpTplName;

        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );
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

void CChisquare2SolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    char tmpChar[256];
    Float64 dtreepathnum;
    store.GetParam( "dtreepathnum", dtreepathnum );
    sprintf(tmpChar, "%.2f", dtreepathnum);

    Float64 redshift;
    Float64 merit;
    Float64 evidence;
    std::string tplName="-1";

    //unused
    Float64 amp;
    Float64 dustCoeff;
    Int32 meiksinIdx;

    if(m_bestRedshiftMethod==0)
    {
        GetBestRedshift( store, redshift, merit, tplName, amp, dustCoeff, meiksinIdx );
        Log.LogInfo( "Chisquare2solve-result: extracted best redshift from chi2 extrema: z=%f", redshift);
    }else if(m_bestRedshiftMethod==2)
    {
        GetBestRedshiftFromPdf( store, redshift, merit, evidence );
        Log.LogInfo( "Chisquare2solve-result: extracted best redshift from PDF: z=%f", redshift);
        GetBestModel(store, redshift, tplName);
        Log.LogInfo( "Chisquare2solve-result: extracted best model: model=%s", tplName.c_str());
    }else{
        Log.LogError( "Chisquare2solve-result: can't parse best redshift estimation method");
    }


    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "Chisquare2Solve_" << tmpChar << "\t"
                << "-1" << "\t" //deltaz
                << m_ReliabilityLabel << "\t"
                << "-1" << "\t"
                << "-1" << "\t"
                << "-1" << "\t"
                << "-1" << "\t"
                << m_TypeLabel << std::endl; //reliability label

}

Bool CChisquare2SolveResult::GetBestRedshift( const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName, Float64& amplitude, Float64& dustCoeff, Int32& meiksinIdx ) const
{
    std::string scopeStr;
    if(m_type == nType_raw){
        scopeStr = "chisquare";
    }else if(m_type == nType_all){
        scopeStr = "chisquare";
    }else if(m_type == nType_noContinuum){
        scopeStr = "chisquare_nocontinuum";
    }else if(m_type == nType_continuumOnly){
        scopeStr = "chisquare_continuum";
    }
    std::string scope = store.GetScope( *this ) + "chisquare2solve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());


    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName = "-1";
    Float64 tmpAmp = 0.0;
    Float64 tmpDustCoeff = 0.0;
    Int32 tmpMeiksinIdx = 0;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->ChiSquare[i] < tmpMerit && meritResult->Status[i] == COperator::nStatus_OK )
            {
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpAmp = meritResult->FitAmplitude[i];
                tmpDustCoeff = meritResult->FitDustCoeff[i];
                tmpMeiksinIdx = meritResult->FitMeiksinIdx[i];
                tmpTplName = (*it).first;
            }
        }
    }


    if( tmpMerit < DBL_MAX )
    {
        redshift = tmpRedshift;
        merit = tmpMerit;
        tplName = tmpTplName;
        amplitude = tmpAmp;
        dustCoeff = tmpDustCoeff;
        meiksinIdx = tmpMeiksinIdx;
        return true;
    }

    return false;

}



/**
 * \brief Searches the best_z = argmax(pdf)
 * output: redshift = argmax(pdf)
 * output: merit = chi2(redshift)
 *
 **/
Bool CChisquare2SolveResult::GetBestRedshiftFromPdf( const CDataStore& store, Float64& redshift, Float64& merit, Float64& evidence ) const
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
        Log.LogError( "GetBestRedshiftFromPdf: no pdf results retrieved from scope: %s", scope_res.c_str());
        return false;
    }


    Float64 tmpProbaLog = -DBL_MAX;
    Float64 tmpRedshift = 0.0;

    for( Int32 i=0; i<logzpdf1d->Redshifts.size(); i++ )
    {
        Float64 probaLog = logzpdf1d->valProbaLog[i];
        if(probaLog>tmpProbaLog)
        {
            tmpProbaLog = probaLog;
            tmpRedshift = logzpdf1d->Redshifts[i];
        }
    }


    redshift = tmpRedshift;
    merit = tmpProbaLog;
    evidence = logzpdf1d->valEvidenceLog;
    return true;
}

Int32 CChisquare2SolveResult::GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
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
 * @brief CChisquare2SolveResult::GetBestModel
 * Find the best tpl-name corresponding to the input arg. z and to the minimum Chisquare
 * @param store
 * @param z
 * @param tplName
 * @return
 */
Int32 CChisquare2SolveResult::GetBestModel(const CDataStore& store, Float64 z, std::string& tplName) const
{
    tplName = "-1";
    Bool foundRedshiftAtLeastOnce = false;

    std::string scopeStr;
    if(m_type == nType_raw){
        scopeStr = "chisquare";
    }else if(m_type == nType_all){
        scopeStr = "chisquare";
    }else if(m_type == nType_noContinuum){
        scopeStr = "chisquare_nocontinuum";
    }else if(m_type == nType_continuumOnly){
        scopeStr = "chisquare_continuum";
    }
    std::string scope = store.GetScope( *this ) + "chisquare2solve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    Float64 tmpMerit = DBL_MAX ;
    std::string tmpTplName = "-1";
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>( (*it).second );
        if(meritResult->ChiSquare.size()<1){
            continue;
        }

        Int32 idx=-1;
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->Status[i] == COperator::nStatus_OK )
            {
                Float64 tmpRedshift = meritResult->Redshifts[i];
                if(tmpRedshift == z){
                    idx = i;
                    foundRedshiftAtLeastOnce = true;
                    break;
                }
            }
        }
        if(idx>-1 && meritResult->ChiSquare[idx]<tmpMerit){
            tmpMerit = meritResult->ChiSquare[idx];
            tmpTplName = (*it).first;
        }

    }

    if(foundRedshiftAtLeastOnce){
        tplName = tmpTplName;
    }else{
        return -1;
    }
    return 1;
}


Bool CChisquare2SolveResult::GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates, Int32 n_candidates) const
{
    Log.LogDebug( "CChisquare2SolveResult:GetRedshiftCandidates" );
    redshiftcandidates.clear();

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
        Log.LogError( "GetRedshiftCandidates: no pdf results retrieved from scope: %s", scope_res.c_str());
        return -1;
    }

    // Find extrema
    TPointList extremumList;
    CExtremum extremum(TFloat64Range(0.0, DBL_MAX), n_candidates);
    extremum.Find( logzpdf1d->Redshifts, logzpdf1d->valProbaLog, extremumList );
    for(Int32 k=0; k<extremumList.size(); k++)
    {
        redshiftcandidates.push_back(extremumList[k].X);
    }

    return true;
}

#include <RedshiftLibrary/method/templatefittingresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>
#include <stdio.h>
#include <float.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/extremum/extremum.h>

using namespace NSEpic;

CTemplateFittingSolveResult::CTemplateFittingSolveResult(const EType type, const std::string scope):
    m_type(type),
    m_scope(scope)
{  
    m_name = m_scope2name[m_scope];
}

void CTemplateFittingSolveResult::preSave(const CDataStore& store)
{
    if(m_bestRedshiftMethod==0)//bestXi2
        GetBestRedshift(store);
    if(m_bestRedshiftMethod==2)//bestPDF
    {
        GetBestRedshiftFromPdf(store);
        Log.LogInfo( "%s-result: extracting best redshift from PDF: z=%f", m_name.c_str(), m_redshift);
        Int32 b = GetBestModel(store, m_redshift);
        if(b==-1){
            Log.LogError(" CTemplateFittingSolveResult::preSave: Couldn't find index of %f", m_redshift);
            throw runtime_error("CTemplateFittingSolveResult::preSave: Couldn't find redshift index. Aborting!");
        }
        Log.LogInfo( "%s-result: extracted best model: model=%s", m_name.c_str(), m_tplName.c_str());

    }
    if(m_tplName=="")
        m_tplName = "Undefined";
}

void CTemplateFittingSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{

  
    stream <<  "#Redshifts\tMerit\tTemplate\tAmplitude\tAmplitudeError\tdustcoeff\tmeiksinidx"<< std::endl;

    stream << m_redshift << "\t"
                << m_merit << "\t"
                << m_tplName << "\t"
                << m_amplitude << "\t"
                << m_amplitudeError << "\t"
                << std::setprecision(4) << m_dustCoeff << "\t"
                << m_meiksinIdx << std::endl;

    stream <<  "#Redshifts\tprobaLog\tevidenceLog\tModel"<< std::endl;
    if(m_bestRedshiftMethod==2)
    {
        stream << m_redshift << "\t"
               << m_merit << "\t"
               << m_evidence << "\t"
               << m_tplName << std::endl;
    }else{
        stream <<  "-1\t-1\t-1"<< std::endl;
    }

    stream << std::endl;
    stream << std::endl;
    std::string detailStr;
    GetBestRedshiftPerTemplateString( store, detailStr);

    stream << detailStr.c_str();
}

Bool CTemplateFittingSolveResult::GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const
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

    std::string scope = m_scope + "." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        Float64 tmpMerit = DBL_MAX ;
        Float64 tmpRedshift = 0.0;
        std::string tmpTplName;

        auto meritResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( (*it).second );
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
//called for saving redshift.csv for Xi2
void CTemplateFittingSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    char tmpChar[256];
    Float64 dtreepathnum;
    store.GetParam( "dtreepathnum", dtreepathnum );
    sprintf(tmpChar, "%.2f", dtreepathnum);

    //preSave(store);

    stream << store.GetSpectrumName() << "\t"
           << store.GetProcessingID() << "\t"
                << m_redshift << "\t"
                << m_merit << "\t"
                << m_tplName << "\t"
                << m_name + "_" << tmpChar << "\t"
                << "-1" << "\t" //deltaz
                << m_ReliabilityLabel << "\t"
                << "-1" << "\t"
                << "-1" << "\t"
                << "-1" << "\t"
                << "-1" << "\t"
                << m_TypeLabel << std::endl; //reliability label

}

Bool CTemplateFittingSolveResult::GetBestRedshift( const CDataStore& store) 
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

    std::string scope = m_scope + "." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName = "-1";
    Float64 tmpAmplitude = 0.0;
    Float64 tmpAmplitudeError = 0.0;
    Float64 tmpDustCoeff = 0.0;
    Int32 tmpMeiksinIdx = 0;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( (*it).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritResult->ChiSquare[i] < tmpMerit && meritResult->Status[i] == COperator::nStatus_OK )
            {
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpAmplitude = meritResult->FitAmplitude[i];
                tmpAmplitudeError = meritResult->FitAmplitudeError[i];
                tmpDustCoeff = meritResult->FitDustCoeff[i];
                tmpMeiksinIdx = meritResult->FitMeiksinIdx[i];
                tmpTplName = (*it).first;
            }
        }
    }

    if( tmpMerit < DBL_MAX )
    {
        m_redshift = tmpRedshift;
        m_merit = tmpMerit;
        m_tplName = tmpTplName;
        m_amplitude = tmpAmplitude;
        m_amplitudeError = tmpAmplitudeError;
        m_dustCoeff = tmpDustCoeff;
        m_meiksinIdx = tmpMeiksinIdx;          
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
Bool CTemplateFittingSolveResult::GetBestRedshiftFromPdf( const CDataStore& store ) 
{
    //check if the pdf is stellar/galaxy or else
    std::string outputPdfRelDir  = "zPDF"; //default value set to galaxy zPDF folder
    std::size_t foundstr_stellar = m_scope.find("stellarsolve");
    std::size_t foundstr_qso = m_scope.find("qsosolve");
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

    m_redshift = tmpRedshift;
    m_merit = tmpProbaLog;
    m_evidence = logzpdf1d->valEvidenceLog;

    return true;
}

Int32 CTemplateFittingSolveResult::GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) 
{
    //check if the pdf is stellar/galaxy or else
    std::string outputPdfRelDir  = "zPDF"; //default value set to galaxy zPDF folder
    std::size_t foundstr_stellar = m_scope.find("stellarsolve");
    std::size_t foundstr_qso = m_scope.find("qsosolve");
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

    m_evidence = logzpdf1d->valEvidenceLog;
    evidence = m_evidence; //prob not useful
    return 0;
}

/**
 * @brief CTemplateFittingSolveResult::GetBestModel
 * Find the best tpl-name corresponding to the input arg. z and to the minimum TemplateFitting
 * @param store
 * @param z
 * @param tplName
 * @param MeiksinIdx (optional return)
 * @param DustCoeff (optional return)
 * @return
 */
Int32 CTemplateFittingSolveResult::GetBestModel(const CDataStore& store, Float64 z)
{
    Bool foundRedshiftAtLeastOnce = false;

    std::string scopeStr;
    if(m_type == nType_raw){
        scopeStr = "templatefitting";
    }else if(m_type == nType_all){
        scopeStr = "templatefitting";
    }else if(m_type == nType_noContinuum){
        scopeStr = "templatefitting_nocontinuum";
    }else if(m_type == nType_continuumOnly){
        scopeStr = "templatefitting_continuum";
    }

    std::string scope = m_scope +"." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    Float64 tmpMerit = DBL_MAX ;
    std::string tmpTplName = "-1";
    Int32 tmpMeiksinIdx = -1;
    Float64 tmpDustCoeff = -1.0;
    Float64 tmpAmplitude = 0.0;
    Float64 tmpAmplitudeError = 0.0;
    Int32 idx = -1;
    bool first = true;
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( (*it).second );
        if(meritResult->ChiSquare.size()<1){
            continue;
        }

       //find redshift index once for all and probably re-search for it if template size is not the same
        if(first || (idx>-1 && meritResult->Redshifts[idx] != z))
        {
           //find the corresponding Z
            auto itZ = std::find(meritResult->Redshifts.begin(), meritResult->Redshifts.end(), z);
            idx = std::distance(meritResult->Redshifts.begin(), itZ);
            first = false;
        }
        if(idx == -1)
            continue;
        //verify that the redshift is valid 
        if(meritResult->Status[idx] == COperator::nStatus_OK) {
            //all conditions are satisfied
            foundRedshiftAtLeastOnce = true;
        }
        if(meritResult->ChiSquare[idx]<tmpMerit){
            tmpMerit = meritResult->ChiSquare[idx];
            tmpTplName = (*it).first;
            tmpMeiksinIdx = meritResult->FitMeiksinIdx[idx];
            tmpDustCoeff = meritResult->FitDustCoeff[idx];
            tmpAmplitude = meritResult->FitAmplitude[idx];
            tmpAmplitudeError = meritResult->FitAmplitudeError[idx];
        }

    }

    if(foundRedshiftAtLeastOnce){
        m_merit = tmpMerit;
        m_tplName = tmpTplName;
        m_meiksinIdx = tmpMeiksinIdx;
        m_dustCoeff = tmpDustCoeff;
        m_amplitude = tmpAmplitude;
        m_amplitudeError = tmpAmplitudeError;
    }else{
        return -1;
    }
    return 1;
}

Bool CTemplateFittingSolveResult::GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates, Int32 n_candidates, std::string outputPdfRelDir) const
{
    Log.LogDebug( "C%sSolveResult::GetRedshiftCandidates", m_name.c_str() );
    redshiftcandidates.clear();

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

void CTemplateFittingSolveResult::getData(const std::string& name, Float64& v) const
{
  if (name.compare("snrHa") == 0)  v = -1;
  else if (name.compare("lfHa") == 0)  v = -1;
  else if (name.compare("snrOII") == 0)  v = -1;
  else if (name.compare("lfOII") == 0)  v = -1;
  else if (name.compare("ContinuumIsmCoeff")== 0)  v = m_dustCoeff;
  else throw Exception("Unknown data %s",name.c_str());
}

void CTemplateFittingSolveResult::getData(const std::string& name, std::string& v) const
{
  if (name.compare("TemplateName") == 0) v = m_tplName;
  else throw Exception("Unknown data %s",name.c_str());
} 

void CTemplateFittingSolveResult::getData(const std::string& name, Int32& v) const
{
  if (name.compare("ContinuumIgmIndex") == 0) v = m_meiksinIdx;
  else throw Exception("Unknown data %s",name.c_str());
} 

void CTemplateFittingSolveResult::getCandidateData(const int& rank,const std::string& name, Float64& v) const
{
  if (name.compare("ContinuumIsmCoeff")== 0)  v = m_dustCoeff;
  else if (name.compare("ContinuumAmplitude") == 0) v = m_amplitude;
  else if (name.compare("ContinuumAmplitudeError") == 0) v = m_amplitudeError; 
  else throw Exception("Unknown data %s",name.c_str());
}

void CTemplateFittingSolveResult::getCandidateData(const int& rank,const std::string& name, std::string& v) const
{
  if (name.compare("TemplateName") == 0) v = m_tplName;
  else throw Exception("Unknown data %s",name.c_str());
} 

void CTemplateFittingSolveResult::getCandidateData(const int& rank,const std::string& name, Int32& v) const
{
  if (name.compare("ContinuumIgmIndex") == 0) v = m_meiksinIdx;
  else throw Exception("Unknown data %s",name.c_str());
} 


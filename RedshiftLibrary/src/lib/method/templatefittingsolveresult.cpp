// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/method/templatefittingsolveresult.h"

#include <stdio.h>
#include <float.h>
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/method/solveresult.h"
#include <memory>

using namespace NSEpic;

CTemplateFittingSolveResult::CTemplateFittingSolveResult(const std::string & scope, 
                                                         const TCandidateZ & BestExtremumResult,
                                                         const std::string & opt_pdfcombination,
                                                         Float64 evidence ):
    CPdfSolveResult( BestExtremumResult, opt_pdfcombination, evidence),
    m_scope(scope)
{}

// TODO DV erase this ?
/*
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

}*/


/*
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
        m_EbmvCoeff = tmpDustCoeff;
        m_meiksinIdx = tmpMeiksinIdx;          
        return true;
    }

    return false;

}
*/

/**
 * \brief Searches the best_z = argmax(pdf)
 * output: redshift = argmax(pdf)
 * output: merit = chi2(redshift)
 *
 **/
/*
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
*/
/*
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
}*/

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
/*
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
        m_EbmvCoeff = tmpDustCoeff;
        m_amplitude = tmpAmplitude;
        m_amplitudeError = tmpAmplitudeError;
    }else{
        return -1;
    }
    return 1;
}
*/


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
#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"

#include "RedshiftLibrary/operator/modelspectrumresult.h"

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/parameterstore.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/statistics/zprior.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/operator/templatefittinglog.h"
#include "RedshiftLibrary/operator/templatefittingwithphot.h"

using namespace NSEpic;
using namespace std;

CMethodTemplateFittingSolve::CMethodTemplateFittingSolve(TScopeStack &scope,string objectType):
  CSolve("templatefittingsolve",scope,objectType)
{
}

std::shared_ptr<CSolveResult> CMethodTemplateFittingSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                                                                   std::shared_ptr<COperatorResultStore> resultStore,
                                                                   TScopeStack &scope)
{

  const CSpectrum& spc=*(inputContext->GetSpectrum());
  const CSpectrum& rebinnedSpc=*(inputContext->GetRebinnedSpectrum());
  const CTemplateCatalog& tplCatalog=*(inputContext->GetTemplateCatalog());
  
  m_redshiftSeparation = inputContext->GetParameterStore()->Get<Float64>( "extremaredshiftseparation");//todo: deci

  Bool storeResult = false;
  Float64 overlapThreshold=inputContext->GetParameterStore()->GetScoped<Float64>( "overlapThreshold");
  std::string opt_spcComponent = inputContext->GetParameterStore()->GetScoped<std::string>( "spectrum.component");
  std::string opt_interp = inputContext->GetParameterStore()->GetScoped<std::string>( "interpolation");
  //disable extinction for stars
  const bool opt_extinction = (m_objectType == "star")? false:inputContext->GetParameterStore()->GetScoped<bool>("extinction");
  bool opt_dustFit = inputContext->GetParameterStore()->GetScoped<bool>("dustfit");

  //std::string calibration_dir = inputContext->GetParameterStore()->Get<std::string>("calibrationDir");
  bool fft_processing = inputContext->GetParameterStore()->GetScoped<bool>("fftprocessing");

  bool use_photometry = false;
  if (inputContext->GetParameterStore()->HasScoped<bool>("enablephotometry"))
    use_photometry = inputContext->GetParameterStore()->GetScoped<bool>("enablephotometry");
  
  if (fft_processing && use_photometry)
    throw GlobalException(INTERNAL_ERROR, "CMethodTemplateFittingSolve::compute: fftprocessing not implemented with photometry enabled");

  if(fft_processing)
    {
        m_templateFittingOperator = std::make_shared<COperatorTemplateFittingLog>(spc, rebinnedSpc, m_lambdaRange, m_redshifts);
        tplCatalog.m_logsampling = true; 
    }
  else{
      if (use_photometry){
        m_templateFittingOperator = std::make_shared<COperatorTemplateFitting>(spc, m_lambdaRange, m_redshifts);
      }else{
        m_templateFittingOperator = std::make_shared<COperatorTemplateFittingPhot>(spc, m_lambdaRange, inputContext->GetPhotBandCatalog(), m_redshifts);
      }
        tplCatalog.m_logsampling = false;
  }


        // prepare the unused masks
    std::vector<CMask> maskList;
        //define the redshift search grid
        //        Log.LogInfo("Stellar fitting redshift range = [%.5f, %.5f], step=%.6f", starRedshiftRange.GetBegin(), starRedshiftRange.GetEnd(), starRedshiftStep);
  
    std::string scopeStr = "templatefitting";
    

    EType _type;
    if(opt_spcComponent=="raw"){
       _type = nType_raw;
    }else if(opt_spcComponent=="nocontinuum"){
       _type = nType_noContinuum;
       scopeStr = "templatefitting_nocontinuum";
    }else if(opt_spcComponent=="continuum"){
        _type = nType_continuumOnly;
        scopeStr = "templatefitting_continuum";
    }else if(opt_spcComponent=="all"){
        _type = nType_all;
    }

    m_opt_maxCandidate = inputContext->GetParameterStore()->GetScoped<int>( "extremacount");
    m_opt_pdfcombination=inputContext->GetParameterStore()->GetScoped<std::string>( "pdfcombination");

    //TODO totaly remove this option ?
    m_opt_enableSaveIntermediateTemplateFittingResults = false;

    Log.LogInfo( "Method parameters:");
    Log.LogInfo( "    -overlapThreshold: %.3f", overlapThreshold);
    Log.LogInfo( "    -component: %s", opt_spcComponent.c_str());
    Log.LogInfo( "    -interp: %s", opt_interp.c_str());
    Log.LogInfo( "    -IGM extinction: %s", opt_extinction?"true":"false");
    Log.LogInfo( "    -ISM dust-fit: %s", opt_dustFit?"true":"false");
    Log.LogInfo( "    -pdfcombination: %s", m_opt_pdfcombination.c_str());
    Log.LogInfo( "    -saveintermediateresults: %d", (int)m_opt_enableSaveIntermediateTemplateFittingResults);
    Log.LogInfo( "");

    Log.LogInfo( "Iterating over %d tplCategories", m_categoryList.size());
    if (tplCatalog.GetTemplateCount(m_categoryList[0]) == 0)
      {
	throw GlobalException(BAD_TEMPLATECATALOG,Formatter()<<"Template catalog for category "<< m_categoryList[0] <<" is empty");
      }
    
    for( UInt32 i=0; i<m_categoryList.size(); i++ )
    {
        std::string category = m_categoryList[i];

        Log.LogInfo( "   trying %s (%d templates)", category.c_str(), tplCatalog.GetTemplateCount( category ));
        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            std::shared_ptr<const CTemplate> tpl = tplCatalog.GetTemplate( category, j );

            Solve(resultStore, fft_processing?rebinnedSpc:spc, tpl, overlapThreshold, maskList, _type, opt_interp, opt_extinction, opt_dustFit);

            storeResult = true;
        }
    }

    if( storeResult )
    {
        COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation, 0.0, m_opt_maxCandidate);   

        std::shared_ptr<PdfCandidatesZResult> candidateResult = pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

        // save in resultstore pdf results
        //        std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
        resultStore->StoreScopedGlobalResult( "pdf", pdfz.m_postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF (deprecated comment, must be removed)
        
        // save in resultstore candidates results
        resultStore->StoreScopedGlobalResult( "candidatesresult", candidateResult );


        //for each extrema, get best model by reading from datastore and selecting best fit
        /////////////////////////////////////////////////////////////////////////////////////
        if (fft_processing){ 
            tplCatalog.m_logsampling = false;
        }
        std::shared_ptr<const ExtremaResult> extremaResult =    
                        SaveExtremaResult( resultStore, scopeStr,
                                               candidateResult->m_ranked_candidates,
                                               tplCatalog,
                                               m_categoryList,
                                               overlapThreshold,
                                               opt_interp);

        // store extrema results
        StoreExtremaResults(resultStore, extremaResult);

        std::shared_ptr< CTemplateFittingSolveResult> TemplateFittingSolveResult = 
                        std::make_shared<CTemplateFittingSolveResult>(resultStore->GetCurrentScopeName(),
                                                                      extremaResult->m_ranked_candidates[0].second,
                                                                      m_opt_pdfcombination,
                                                                      pdfz.m_postmargZResult->valEvidenceLog);

        return TemplateFittingSolveResult;
    }

    return NULL;
}


Bool CMethodTemplateFittingSolve::Solve(std::shared_ptr<COperatorResultStore> resultStore,
                                   const CSpectrum& spc,
                                   const std::shared_ptr<const CTemplate> & tpl,
                                   Float64 overlapThreshold,
                                   std::vector<CMask> maskList,
                                   EType spctype,
                                   std::string opt_interp,
                                   bool opt_extinction,
                                   bool opt_dustFitting)
{

    std::string scopeStr = "templatefitting";
    Int32 _ntype = 1;
    CSpectrum::EType _spctype = CSpectrum::nType_raw;
    CSpectrum::EType _spctypetab[3] = {CSpectrum::nType_raw, CSpectrum::nType_noContinuum, CSpectrum::nType_continuumOnly};

    Int32 enable_extinction = 0; //TODO: extinction should be deactivated for nocontinuum anyway ? TBD
    if(opt_extinction)
    {
        enable_extinction = 1;
    }

    Int32 option_dustFitting = -1;
    if(opt_dustFitting)
    {
        option_dustFitting = -10;
    }

    //case: nType_all
    if(spctype == nType_all){
        _ntype = 3;
    }

    const CSpectrum::EType save_spcType = spc.GetType();
    const CSpectrum::EType save_tplType = tpl->GetType();

    for( Int32 i=0; i<_ntype; i++){
        if(spctype == nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = static_cast<CSpectrum::EType>(spctype);
        }
        spc.SetType(_spctype);
        tpl->SetType(_spctype);

        if(_spctype == CSpectrum::nType_continuumOnly){
            // use continuum only
            scopeStr = "templatefitting_continuum";

        }else if(_spctype == CSpectrum::nType_raw){
            // use full spectrum
            scopeStr = "templatefitting";

        }else if(_spctype == CSpectrum::nType_noContinuum){
            // use spectrum without continuum
            scopeStr = "templatefitting_nocontinuum";
            //
            option_dustFitting = -1;
        }

        // Compute merit function
        auto  templateFittingResult = std::dynamic_pointer_cast<CTemplateFittingResult>( m_templateFittingOperator->Compute( tpl,
                                                                                                                             overlapThreshold, maskList, opt_interp,
                                                                                                                             enable_extinction, option_dustFitting ) );
        
        if( !templateFittingResult )
        {
            //Log.LogError( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            resultStore->StoreScopedPerTemplateResult( tpl, scopeStr.c_str(), templateFittingResult );

            //Save intermediate templatefitting results
            if(m_opt_enableSaveIntermediateTemplateFittingResults && templateFittingResult->ChiSquareIntermediate.size()>0 && templateFittingResult->ChiSquareIntermediate.size()==templateFittingResult->Redshifts.size())
            {
                Int32 nISM = templateFittingResult->ChiSquareIntermediate[0].size();
                if(templateFittingResult->ChiSquareIntermediate[0].size()>0)
                {
                    Int32 nIGM = templateFittingResult->ChiSquareIntermediate[0][0].size();

                    for(Int32 kism=0; kism<nISM; kism++)
                    {
                        for(Int32 kigm=0; kigm<nIGM; kigm++)
                        {
                            std::shared_ptr<CTemplateFittingResult> result_chisquare_intermediate = std::shared_ptr<CTemplateFittingResult>( new CTemplateFittingResult() );
                            result_chisquare_intermediate->Init( templateFittingResult->Redshifts.size(), 0, 0);
                            for(Int32 kz=0; kz<templateFittingResult->Redshifts.size(); kz++)
                            {
                                result_chisquare_intermediate->Redshifts[kz] = templateFittingResult->Redshifts[kz];
                                result_chisquare_intermediate->ChiSquare[kz] = templateFittingResult->ChiSquareIntermediate[kz][kism][kigm];
                            }

                            std::string resname = (boost::format("%s_intermediate_ism%d_igm%d") % scopeStr.c_str() % kism % kigm).str();
                            resultStore->StoreScopedPerTemplateResult( tpl, resname.c_str(), result_chisquare_intermediate );
                        }
                    }
                }
            }
        }
    }

    spc.SetType(save_spcType);
    tpl->SetType(save_tplType);

    return true;
}

ChisquareArray CMethodTemplateFittingSolve::BuildChisquareArray(std::shared_ptr<const COperatorResultStore> store, const std::string & scopeStr) const
{
    ChisquareArray chisquarearray;

    Log.LogDetail("templatefittingsolve: building chisquare array");
    std::string scope = store->GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    templatefittingsolve: using results in scope: %s", scope.c_str());

    TOperatorResultMap meritResults = store->GetPerTemplateResult(scope.c_str());

    chisquarearray.cstLog = -1;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( (*it).second );
        Int32 nISM = -1;
        Int32 nIGM = -1;
        if(meritResult->ChiSquareIntermediate.size()>0)
        {
            nISM = meritResult->ChiSquareIntermediate[0].size();
            if(meritResult->ChiSquareIntermediate[0].size()>0)
            {
                nIGM = meritResult->ChiSquareIntermediate[0][0].size();
            }
        }
        if(chisquarearray.cstLog==-1)
        {
            chisquarearray.cstLog = meritResult->CstLog;
            Log.LogInfo("templatefittingsolve: using cstLog = %f", chisquarearray.cstLog);
        }else if ( chisquarearray.cstLog != meritResult->CstLog)
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"templatefittingsolve: Found different cstLog values in results... val-1="<<chisquarearray.cstLog <<" != val-2="<< meritResult->CstLog);
        }
        if(chisquarearray.redshifts.size()==0)
        {
            chisquarearray.redshifts = meritResult->Redshifts;
        }

        //check chi2 results status for this template
        {
            Bool foundBadStatus = 0;
            for ( UInt32 kz=0; kz<meritResult->Redshifts.size(); kz++)
            {
                if(meritResult->Status[kz]!=COperator::nStatus_OK)
                {
                    foundBadStatus = 1;
                    break;
                }
            }
            if(foundBadStatus)
            {
	      throw GlobalException(INTERNAL_ERROR,Formatter()<<"templatefittingsolve: Found bad status result... for tpl="<< (*it).first.c_str());
            }
        }

        CZPrior zpriorhelper;
        for(Int32 kism=0; kism<nISM; kism++)
        {
            for(Int32 kigm=0; kigm<nIGM; kigm++)
            {
                chisquarearray.zpriors.push_back(zpriorhelper.GetConstantLogZPrior(meritResult->Redshifts.size()));

                //correct chi2 for ampl. marg. if necessary: todo add switch, currently deactivated
                chisquarearray.chisquares.emplace_back(meritResult->ChiSquareIntermediate.size(), DBL_MAX);
                TFloat64List & logLikelihoodCorrected = chisquarearray.chisquares.back();
                for ( UInt32 kz=0; kz<meritResult->Redshifts.size(); kz++)
                {
                    logLikelihoodCorrected[kz] = meritResult->ChiSquareIntermediate[kz][kism][kigm];// + resultXXX->ScaleMargCorrectionTplshapes[][]?;
                }
                Log.LogDetail("    templatefittingsolve: Pdfz combine - prepared merit  #%d for model : %s", chisquarearray.chisquares.size()-1, ((*it).first).c_str());
            }
        }
    }

    return chisquarearray;
}


std::shared_ptr<const ExtremaResult> 
CMethodTemplateFittingSolve::SaveExtremaResult(std::shared_ptr<const COperatorResultStore> store,
                                               const std::string & scopeStr,
                                               const TCandidateZbyRank & ranked_zCandidates,
                                               const CTemplateCatalog& tplCatalog,
                                               const TStringList& tplCategoryList,
                                               Float64 overlapThreshold,
                                               std::string opt_interp)
{

    Log.LogDetail("CMethodTemplateFittingSolve::SaveExtremaResult: building chisquare array");
    std::string scope = store->GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    templatefittingsolve: using results in scope: %s", scope.c_str());

    TOperatorResultMap results = store->GetPerTemplateResult(scope.c_str());

    //Bool foundRedshiftAtLeastOnce = false;

    Int32 extremumCount = ranked_zCandidates.size();
    
    auto firstResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( (*results.begin()).second );
    const TFloat64List & redshifts = firstResult->Redshifts;

    //check all results  status
    for ( auto & r : results)
    {
        auto TplFitResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( r.second );
        if(TplFitResult->ChiSquare.size() != redshifts.size()){
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"CMethodTemplateFittingSolve::SaveExtremaResult, templatefitting results (for tpl="<< r.first.c_str()<<") has wrong size");
        }

        Bool foundBadStatus = false;
        Bool foundBadRedshift = false;
        for ( UInt32 kz=0; kz<TplFitResult->Redshifts.size(); kz++)
        {
            if(TplFitResult->Status[kz]!=COperator::nStatus_OK)
            {
                foundBadStatus = true;
                break;
            }
            if(TplFitResult->Redshifts[kz]!=redshifts[kz])
            {
                foundBadRedshift = true;
                break;
            }
        }
        if(foundBadStatus)
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"CMethodTemplateFittingSolve::SaveExtremaResult: Found bad status result... for tpl="<< r.first.c_str());
        }
        if(foundBadRedshift)
        {
	  throw GlobalException(INTERNAL_ERROR,Formatter()<<"CMethodTemplateFittingSolve::SaveExtremaResult: redshift vector is not the same for tpl="<< r.first.c_str());
        }

    }

    std::shared_ptr<ExtremaResult> extremaResult = make_shared<ExtremaResult>(ranked_zCandidates);


    for (Int32 i = 0; i < extremumCount; i++)
    {
        //std::string Id = ranked_zCandidates[i].first;
        Float64 z = ranked_zCandidates[i].second.Redshift;

        //find the corresponding Z
        auto itZ = std::find(redshifts.begin(), redshifts.end(), z);
        const Int32 idx = std::distance(redshifts.begin(), itZ);

        // find the min chisquare at corresponding redshift
        Float64 ChiSquare = DBL_MAX ;
        std::string tplName = "";
        for( auto & r : results)
        {
            auto TplFitResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( r.second );

            if(TplFitResult->ChiSquare[idx] < ChiSquare) {
                ChiSquare = TplFitResult->ChiSquare[idx];
                tplName = r.first;
            };
        }

        // Fill extrema Result
        auto TplFitResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( results[tplName] );
        extremaResult->m_ranked_candidates[i].second.FittedTplMerit = ChiSquare;
        extremaResult->m_ranked_candidates[i].second.FittedTplName= tplName;
        extremaResult->m_ranked_candidates[i].second.FittedTplMeiksinIdx= TplFitResult->FitMeiksinIdx[idx];
        extremaResult->m_ranked_candidates[i].second.FittedTplEbmvCoeff= TplFitResult->FitEbmvCoeff[idx];
        extremaResult->m_ranked_candidates[i].second.FittedTplAmplitude= TplFitResult->FitAmplitude[idx];
        extremaResult->m_ranked_candidates[i].second.FittedTplAmplitudeError= TplFitResult->FitAmplitudeError[idx];
        extremaResult->m_ranked_candidates[i].second.FittedTplDtm= TplFitResult->FitDtM[idx];
        extremaResult->m_ranked_candidates[i].second.FittedTplMtm= TplFitResult->FitMtM[idx];
        extremaResult->m_ranked_candidates[i].second.FittedTplLogPrior= TplFitResult->LogPrior[idx];

        Float64 FitSNR = NAN;
        if (TplFitResult->FitMtM[idx] != 0.)
            FitSNR = abs(TplFitResult->FitDtM[idx])/sqrt(TplFitResult->FitMtM[idx]); // = |amplitude|/amplitudeError
        extremaResult->m_ranked_candidates[i].second.FittedTplSNR= FitSNR;
        //make sure tpl is non-rebinned
        Bool currentSampling = tplCatalog.m_logsampling;
        tplCatalog.m_logsampling=false;
        std::shared_ptr<const CTemplate> tpl = tplCatalog.GetTemplateByName(tplCategoryList, tplName);
        std::shared_ptr<CModelSpectrumResult> spcmodelPtr = 
                            m_templateFittingOperator->ComputeSpectrumModel(tpl, 
                                                                            z,
                                                                            TplFitResult->FitEbmvCoeff[idx],
                                                                            TplFitResult->FitMeiksinIdx[idx],
                                                                            TplFitResult->FitAmplitude[idx],
                                                                            opt_interp, 
                                                                            overlapThreshold);
        
        if(spcmodelPtr==nullptr)
            throw GlobalException(INTERNAL_ERROR,"CMethodTemplateFittingSolve::SaveExtremaResult: Couldnt compute spectrum model");
        tplCatalog.m_logsampling = currentSampling;                                                
        extremaResult->m_savedModelSpectrumResults[i] = std::move(spcmodelPtr);

     
    }

    return extremaResult;
}


void CMethodTemplateFittingSolve::StoreExtremaResults( std::shared_ptr<COperatorResultStore> resultStore, 
                                                       std::shared_ptr<const ExtremaResult> & extremaResult) const
{
  resultStore->StoreScopedGlobalResult("extrema_results",extremaResult);
  Log.LogInfo("CMethodTemplateFittingSolve::StoreExtremaResults: Templatefitting, saving extrema results");
   
  return;
}

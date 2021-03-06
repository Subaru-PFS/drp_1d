#include "RedshiftLibrary/method/templatefittingsolve.h"
#include "RedshiftLibrary/method/templatefittingsolveresult.h"

#include "RedshiftLibrary/operator/modelcontinuumfittingresult.h"
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

using namespace NSEpic;
using namespace std;

CMethodTemplateFittingSolve::CMethodTemplateFittingSolve(TScopeStack &scope,string objectType):
  CSolve("templatefittingsolve",scope,objectType)
{
}

void CMethodTemplateFittingSolve::GetRedshiftSampling(std::shared_ptr<const CInputContext> inputContext, TFloat64Range& redshiftRange, Float64& redshiftStep)
{
    if(m_objectType == "galaxy" && inputContext->m_use_LogLambaSpectrum)
    {
       redshiftRange = inputContext->m_redshiftRangeFFT;
       redshiftStep = inputContext->m_redshiftStepFFT;
       if(m_redshiftSampling=="lin"){
            m_redshiftSampling = "log";
            Log.LogWarning("m_redshift sampling value is forced to log since FFTprocessing is used");
       }
    }else
    {        
        //default is to read from the scoped paramStore
        CSolve::GetRedshiftSampling(inputContext, redshiftRange, redshiftStep );
    }
}

std::shared_ptr<CSolveResult> CMethodTemplateFittingSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                                                                   std::shared_ptr<COperatorResultStore> resultStore,
                                                                   TScopeStack &scope)
{

  const CSpectrum& spc=*(inputContext->GetSpectrum().get());
  const CSpectrum& rebinnedSpc=*(inputContext->GetRebinnedSpectrum().get());
  const CTemplateCatalog& tplCatalog=*(inputContext->GetTemplateCatalog().get());
  
  m_redshiftSeparation = inputContext->GetParameterStore()->Get<Float64>( "extremaredshiftseparation");//todo: deci

  Bool storeResult = false;
  Float64 overlapThreshold=inputContext->GetParameterStore()->GetScoped<Float64>( "overlapThreshold");
  std::string opt_spcComponent = inputContext->GetParameterStore()->GetScoped<std::string>( "spectrum.component");
  std::string opt_interp = inputContext->GetParameterStore()->GetScoped<std::string>( "interpolation");
  //disable extinction for stars
  const std::string opt_extinction = (m_objectType == "star")? "no":inputContext->GetParameterStore()->GetScoped<std::string>("extinction");
  std::string opt_dustFit = inputContext->GetParameterStore()->GetScoped<std::string>("dustfit");

  //std::string calibration_dir = inputContext->GetParameterStore()->Get<std::string>("calibrationDir");
  bool fft_processing = inputContext->GetParameterStore()->GetScoped<std::string>("fftprocessing") == "yes";
  
  if(fft_processing)
    {
      //  std::string opt_interp="lin";//this is important for creating the spectrum model//moved from processflow
      m_templateFittingOperator = std::shared_ptr<COperatorTemplateFittingBase>(new COperatorTemplateFittingLog());
      tplCatalog.m_logsampling = true; 
    }
  else{
       m_templateFittingOperator = std::shared_ptr<COperatorTemplateFittingBase>(new COperatorTemplateFitting());
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
    /*
    m_opt_saveintermediateresults = inputContext->GetParameterStore()->GetScoped<std::string>( "saveintermediateresults");
    if(m_opt_saveintermediateresults=="yes")
    {
        m_opt_enableSaveIntermediateTemplateFittingResults = true;
    }else{
        m_opt_enableSaveIntermediateTemplateFittingResults = false;
    }
    */
    //TODO totaly remove this option ?
    m_opt_enableSaveIntermediateTemplateFittingResults = false;

    Log.LogInfo( "Method parameters:");
    Log.LogInfo( "    -overlapThreshold: %.3f", overlapThreshold);
    Log.LogInfo( "    -component: %s", opt_spcComponent.c_str());
    Log.LogInfo( "    -interp: %s", opt_interp.c_str());
    Log.LogInfo( "    -IGM extinction: %s", opt_extinction.c_str());
    Log.LogInfo( "    -ISM dust-fit: %s", opt_dustFit.c_str());
    Log.LogInfo( "    -pdfcombination: %s", m_opt_pdfcombination.c_str());
    Log.LogInfo( "    -saveintermediateresults: %d", (int)m_opt_enableSaveIntermediateTemplateFittingResults);
    Log.LogInfo( "");

    Log.LogInfo( "Iterating over %d tplCategories", m_categoryList.size());
    for( UInt32 i=0; i<m_categoryList.size(); i++ )
    {
        std::string category = m_categoryList[i];

        Log.LogInfo( "   trying %s (%d templates)", category.c_str(), tplCatalog.GetTemplateCount( category ));
        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            std::shared_ptr<const CTemplate> tpl = tplCatalog.GetTemplate( category, j );

            Solve(resultStore, fft_processing?rebinnedSpc:spc, *tpl, m_lambdaRange, m_redshifts, overlapThreshold, maskList, _type, opt_interp, opt_extinction, opt_dustFit);

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
        //using clamped lambdaRange:
        TFloat64Range clampedLbdaRange;
        spc.GetSpectralAxis().ClampLambdaRange( m_lambdaRange, clampedLbdaRange );
        std::shared_ptr<const ExtremaResult> extremaResult = 
                        SaveExtremaResult( resultStore, scopeStr,
                                               candidateResult->m_ranked_candidates,
                                               spc,
                                               tplCatalog,
                                               m_categoryList,
                                               clampedLbdaRange,
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
                                   const CTemplate& tpl,
                                   const TFloat64Range& lambdaRange,
                                   const TFloat64List& redshifts,
                                   Float64 overlapThreshold,
                                   std::vector<CMask> maskList,
                                   EType spctype,
                                   std::string opt_interp,
                                   std::string opt_extinction,
                                   std::string opt_dustFitting)
{

    std::string scopeStr = "templatefitting";
    Int32 _ntype = 1;
    CSpectrum::EType _spctype = CSpectrum::nType_raw;
    CSpectrum::EType _spctypetab[3] = {CSpectrum::nType_raw, CSpectrum::nType_noContinuum, CSpectrum::nType_continuumOnly};

    Int32 enable_extinction = 0; //TODO: extinction should be deactivated for nocontinuum anyway ? TBD
    if(opt_extinction=="yes")
    {
        enable_extinction = 1;
    }

    Int32 option_dustFitting = -1;
    if(opt_dustFitting=="yes")
    {
        option_dustFitting = -10;
    }

    //case: nType_all
    if(spctype == nType_all){
        _ntype = 3;
    }

    const CSpectrum::EType save_spcType = spc.GetType();
    const CSpectrum::EType save_tplType = tpl.GetType();

    for( Int32 i=0; i<_ntype; i++){
        if(spctype == nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = static_cast<CSpectrum::EType>(spctype);
        }
        spc.SetType(_spctype);
        tpl.SetType(_spctype);

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
        //CRef<CTemplateFittingResult>  templateFittingResult = (CTemplateFittingResult*)templatefitting.ExportChi2versusAZ( _spc, _tpl, lambdaRange, redshifts, overlapThreshold );
        auto  templateFittingResult = std::dynamic_pointer_cast<CTemplateFittingResult>( m_templateFittingOperator->Compute( spc, tpl,lambdaRange,redshifts,
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
    tpl.SetType(save_tplType);

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
            Log.LogError("templatefittingsolve: Found different cstLog values in results... val-1=%f != val-2=%f", chisquarearray.cstLog, meritResult->CstLog);
            throw runtime_error("templatefittingsolve: Found different cstLog values in results");
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
                Log.LogError("templatefittingsolve: Found bad status result... for tpl=%s", (*it).first.c_str());
                throw runtime_error("templatefittingsolve: Found bad status result");
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
                                               const CSpectrum& spc,
                                               const CTemplateCatalog& tplCatalog,
                                               const TStringList& tplCategoryList,
                                               const TFloat64Range& lambdaRange,
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
            Log.LogError("CMethodTemplateFittingSolve::SaveExtremaResult, templatefitting results (for tpl=%s) has wrong size", r.first.c_str());
            throw runtime_error("CMethodTemplateFittingSolve::SaveExtremaResult, one templatefitting results has wrong size");
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
            Log.LogError("CMethodTemplateFittingSolve::SaveExtremaResult: Found bad status result... for tpl=%s", r.first.c_str());
            throw runtime_error("CMethodTemplateFittingSolve::SaveExtremaResult: Found bad status result");
        }
        if(foundBadRedshift)
        {
            Log.LogError("CMethodTemplateFittingSolve::SaveExtremaResult: redshift vector is not the same for tpl=%s", r.first.c_str());
            throw runtime_error("CMethodTemplateFittingSolve::SaveExtremaResult: Found different redshift vector");
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
        std::shared_ptr<CModelSpectrumResult> spcmodelPtr; 
        m_templateFittingOperator->ComputeSpectrumModel(spc, *tpl, 
                                                        z,
                                                        TplFitResult->FitEbmvCoeff[idx],
                                                        TplFitResult->FitMeiksinIdx[idx],
                                                        TplFitResult->FitAmplitude[idx],
                                                        opt_interp, lambdaRange, 
                                                        overlapThreshold, spcmodelPtr);
        tplCatalog.m_logsampling = currentSampling;                                                
        extremaResult->m_savedModelSpectrumResults[i] = std::move(spcmodelPtr);

        extremaResult->m_savedModelContinuumFittingResults[i] = 
                        std::make_shared<CModelContinuumFittingResult>( z,
                                                                        tplName,
                                                                        ChiSquare,
                                                                        TplFitResult->FitAmplitude[idx],
                                                                        TplFitResult->FitAmplitudeError[idx],
                                                                        TplFitResult->FitEbmvCoeff[idx],
                                                                        TplFitResult->FitMeiksinIdx[idx],
                                                                        FitSNR );          
    }

    return extremaResult;
}


void CMethodTemplateFittingSolve::StoreExtremaResults( std::shared_ptr<COperatorResultStore> resultStore, 
                                                       std::shared_ptr<const ExtremaResult> & extremaResult) const
{
  resultStore->StoreScopedGlobalResult("extrema_results",extremaResult);
  Log.LogInfo("Linemodel, saving extrema results");
   
  return;
}

#include <RedshiftLibrary/method/tplcombinationsolve.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/statistics/zprior.h>

#include <cfloat>

using namespace NSEpic;
using namespace std;


CMethodTplcombinationSolve::CMethodTplcombinationSolve(TScopeStack &scope,string objectType):
  CSolve("tplcombinationsolve",scope,objectType)
{
}



std::shared_ptr<CSolveResult> CMethodTplcombinationSolve::compute(std::shared_ptr<const CInputContext> inputContext,
                                                                                 std::shared_ptr<COperatorResultStore> resultStore,
                                                                                 TScopeStack &scope)


{

  const CSpectrum& spc=*(inputContext->GetSpectrum().get());
  const CTemplateCatalog& tplCatalog=*(inputContext->GetTemplateCatalog().get());
  const CRayCatalog& restraycatalog=*(inputContext->GetRayCatalog().get());

    Bool storeResult = false;
    m_redshiftSeparation = inputContext->GetParameterStore()->Get<Float64>("extremaredshiftseparation");
    m_opt_maxCandidate = inputContext->GetParameterStore()->GetScoped<int>( "extremacount");
    m_opt_pdfcombination = inputContext->GetParameterStore()->GetScoped<std::string>( "pdfcombination");
    std::string opt_interp = inputContext->GetParameterStore()->GetScoped<std::string>( "interpolation");
  std::string opt_dustFit = inputContext->GetParameterStore()->GetScoped<std::string>("dustfit");
  Float64 overlapThreshold=inputContext->GetParameterStore()->GetScoped<Float64>( "overlapThreshold");
  std::string opt_extinction = inputContext->GetParameterStore()->GetScoped<std::string>("extinction");
 
    m_opt_saveintermediateresults = inputContext->GetParameterStore()->GetScoped<std::string>( "saveintermediateresults");
    std::string opt_spcComponent = inputContext->GetParameterStore()->GetScoped<std::string>( "spectrum.component");

    std::vector<CMask> maskList;    

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

    if(m_opt_saveintermediateresults=="yes")
    {
        m_opt_enableSaveIntermediateChisquareResults = true;
    }else{
        m_opt_enableSaveIntermediateChisquareResults = false;
    }

    //for now interp must be 'lin'. pfg not availbale for now...
    if(opt_interp!="lin")
    {
        Log.LogError("Tplcombinationsolve: interp. parameter must be 'lin'");
        throw runtime_error("Tplcombinationsolve: interpolation parameter must be lin");
    }

    Log.LogInfo( "Method parameters:");
    Log.LogInfo( "    -interpolation: %s", opt_interp.c_str());
    Log.LogInfo( "    -overlapThreshold: %.3f", overlapThreshold);
    Log.LogInfo( "    -component: %s", opt_spcComponent.c_str());
    Log.LogInfo( "    -IGM extinction: %s", opt_extinction.c_str());
    Log.LogInfo( "    -ISM dust-fit: %s", opt_dustFit.c_str());
    //Log.LogInfo( "    -pdfcombination: %s", m_opt_pdfcombination.c_str());
    Log.LogInfo( "    -saveintermediateresults: %d", (int)m_opt_enableSaveIntermediateChisquareResults);
    Log.LogInfo( "");



    
    Solve(resultStore,
          spc,
          tplCatalog,
          m_categoryList,
          m_lambdaRange,
          m_redshifts,
          overlapThreshold,
          maskList,
          _type,
          opt_interp,
          opt_extinction,
          opt_dustFit);



    COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation, 0.0, m_opt_maxCandidate);

    std::shared_ptr<CPdfCandidateszResult> candidateResult = pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

    // save in resultstore pdf results

    resultStore->StoreGlobalResult( "pdf", pdfz.m_postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF

    // save in resultstore candidates results
    resultStore->StoreScopedGlobalResult("candidatesresult", candidateResult );


    //for each candidate, get best model by reading from datastore and selecting best fit
    /////////////////////////////////////////////////////////////////////////////////////

    
    std::shared_ptr<const CExtremaResult> ExtremaResult = std::make_shared<CExtremaResult>();
    //TODO this is taken from templatefittingsolve, make it work here ?
    // common base class for templatefitting and tplcombination 
      /*                    SaveExtremaResult( resultStore, scopeStr,
                                               candidateResult->m_ranked_candidates,
                                               spc,
                                               tplCatalog,
                                               m_categoryList,
                                               m_lambdaRange,
                                               overlapThreshold,
                                               opt_interp,
                                               opt_extinction );
    // store extrema results
        StoreExtremaResults(resultStore, ExtremaResult);
    */


    
    std::shared_ptr< CTemplateFittingSolveResult> solveResult =
      std::make_shared< CTemplateFittingSolveResult>(resultStore->GetCurrentScopeName(),
                                                     ExtremaResult,
                                                     m_opt_pdfcombination,
                                                     pdfz.m_postmargZResult->valEvidenceLog);

    // TBD

    return solveResult;

}

Bool CMethodTplcombinationSolve::Solve(std::shared_ptr<COperatorResultStore> resultStore,
                                       const CSpectrum& spc,
                                       const CTemplateCatalog& tplCatalog,
                                       const TStringList& tplCategoryList,
                                       const TFloat64Range& lambdaRange,
                                       const TFloat64List& redshifts,
                                       Float64 overlapThreshold,
                                       std::vector<CMask> maskList,
                                       EType spctype,
                                       std::string opt_interp,
                                       std::string opt_extinction,
                                       std::string opt_dustFitting)
{
    std::string scopeStr = "tplcombination";
    Int32 _ntype = 1;
    CSpectrum::EType _spctype = CSpectrum::nType_raw;
    CSpectrum::EType _spctypetab[3] = {CSpectrum::nType_raw, CSpectrum::nType_noContinuum, CSpectrum::nType_continuumOnly};

    Int32 enable_extinction = 0; //TODO: extinction should be deactivated for nocontinuum anyway ? TBD
    if(opt_extinction=="yes")
    {
        enable_extinction = 1;
    }
    Int32 enable_dustFitting = 0;
    if(opt_dustFitting=="yes")
    {
        enable_dustFitting = 1;
    }

    //prepare the list of components/templates
    std::vector<CTemplate> tplList;
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            std::shared_ptr<const CTemplate> tpl = tplCatalog.GetTemplate( category, j );
            tplList.push_back(*tpl);
        }
    }


    //case: nType_all
    if(spctype == nType_all){
        _ntype = 3;
    }

    const CSpectrum::EType save_spcType = spc.GetType();

    std::vector<CSpectrum::EType> save_tplTypes;
    for ( CTemplate tpl: tplList){
        save_tplTypes.push_back(tpl.GetType());
    }

    for( Int32 i=0; i<_ntype; i++){
        if(spctype == nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = static_cast<CSpectrum::EType>(spctype);
        }

        spc.SetType(_spctype);
        for (CTemplate tpl: tplList)
            tpl.SetType(_spctype);

        if(_spctype == CSpectrum::nType_continuumOnly){
            // use continuum only
            scopeStr = "tplcombination_continuum";

        }else if(_spctype == CSpectrum::nType_raw){
            // use full spectrum
            scopeStr = "tplcombination";

        }else if(_spctype == CSpectrum::nType_noContinuum){
            // use spectrum without continuum
            scopeStr = "tplcombination_nocontinuum";
            //
            enable_dustFitting = 0;
        }

        // Compute merit function
        auto  result = std::dynamic_pointer_cast<CTemplateFittingResult>( m_tplcombinationOperator.Compute( spc, 
                                                                                                            tplList, 
                                                                                                            lambdaRange, 
                                                                                                            redshifts, 
                                                                                                            overlapThreshold, 
                                                                                                            maskList, 
                                                                                                            opt_interp, 
                                                                                                            enable_extinction, 
                                                                                                            enable_dustFitting ) );

        if( !result )
        {
            //Log.LogError( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            Log.LogDetail("tplcombinationsolve: Save tplcombination results");
            resultStore->StoreScopedGlobalResult(scopeStr.c_str(), result );
            // Store spectrum results
            Log.LogDetail("tplcombinationsolve: Save spectrum/model results");
            m_tplcombinationOperator.SaveSpectrumResults(resultStore);
        }
    }

    // restore component types
    spc.SetType(save_spcType);
    auto it = save_tplTypes.begin();
    for (CTemplate tpl: tplList)
    {
        tpl.SetType(*it);
        it++;
    }

    return true;
}

ChisquareArray CMethodTplcombinationSolve::BuildChisquareArray(std::shared_ptr<COperatorResultStore> store,
                                                               const std::string & scopeStr) const
{
    ChisquareArray chisquarearray;

    Log.LogDetail("tplcombinationsolve: build chisquare array");
    std::string scope = store->GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    tplcombinationsolve: using results in scope: %s", scope.c_str());

    auto results = store->GetGlobalResult( scope.c_str() );
    if(results.expired())
    {
        throw runtime_error("tplcombinationsolve: CombinePDF - Unable to retrieve tplcombination results");
    }
    std::shared_ptr<const CTemplateFittingResult> result = std::dynamic_pointer_cast<const CTemplateFittingResult>( results.lock() );

    chisquarearray.cstLog = -1;

    Int32 retPdfz=-1;
    {
        Int32 nISM = -1;
        Int32 nIGM = -1;
        if(result->ChiSquareIntermediate.size()>0)
        {
            nISM = result->ChiSquareIntermediate[0].size();
            if(result->ChiSquareIntermediate[0].size()>0)
            {
                nIGM = result->ChiSquareIntermediate[0][0].size();
            }
        }
        if(chisquarearray.cstLog==-1)
        {
            chisquarearray.cstLog = result->CstLog;
            Log.LogInfo("tplcombinationsolve: using cstLog = %f", chisquarearray.cstLog);
        }else if ( chisquarearray.cstLog != result->CstLog)
        {
            Log.LogError("tplcombinationsolve: Found different cstLog values in results... val-1=%f != val-2=%f", chisquarearray.cstLog, result->CstLog);
            throw runtime_error("tplcombinationsolve: Found different cstLog values in results");
        }
        if(chisquarearray.redshifts.size()==0)
        {
            chisquarearray.redshifts = result->Redshifts;
        }

        //check chi2 results status for this template
        {
            Bool foundBadStatus = 0;
            for ( UInt32 kz=0; kz<result->Redshifts.size(); kz++)
            {
                if(result->Status[kz]!=COperator::nStatus_OK)
                {
                    foundBadStatus = 1;
                    break;
                }
            }
            if(foundBadStatus)
            {
                Log.LogError("tplcombinationsolve: Found bad status result...");
                throw runtime_error("tplcombinationsolve: Found bad status result");
            }
        }

        CZPrior zpriorhelper;
        for(Int32 kism=0; kism<nISM; kism++)
        {
            for(Int32 kigm=0; kigm<nIGM; kigm++)
            {
                chisquarearray.zpriors.push_back(zpriorhelper.GetConstantLogZPrior(result->Redshifts.size()));

                //correct chi2 for ampl. marg. if necessary: todo add switch, currently deactivated
                chisquarearray.chisquares.emplace_back(result->ChiSquareIntermediate.size(), DBL_MAX);
                TFloat64List & logLikelihoodCorrected = chisquarearray.chisquares.back();
                for ( UInt32 kz=0; kz<result->Redshifts.size(); kz++)
                {
                    logLikelihoodCorrected[kz] = result->ChiSquareIntermediate[kz][kism][kigm];// + resultXXX->ScaleMargCorrectionTplshapes[][]?;
                }
                Log.LogDetail("    tplcombinationsolve: Pdfz combine - prepared merit #%d for ism=%d, igm=%d", chisquarearray.chisquares.size()-1, kism, kigm);
            }
        }
    }

    return chisquarearray;
}


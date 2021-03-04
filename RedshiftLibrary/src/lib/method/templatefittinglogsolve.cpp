#include <RedshiftLibrary/method/templatefittinglogsolve.h>
#include <RedshiftLibrary/method/templatefittingsolveresult.h>
#include <RedshiftLibrary/operator/templatefitting.h> //needed to compute model at the end
#include <RedshiftLibrary/operator/modelcontinuumfittingresult.h>
#include <RedshiftLibrary/operator/modelspectrumresult.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/statistics/zprior.h>

//#include <cfloat>

using namespace NSEpic;
using namespace std;


CMethodTemplateFittingLogSolve::CMethodTemplateFittingLogSolve( std::string calibrationPath ):
    m_templateFittingOperator( calibrationPath )
{
}

const std::string CMethodTemplateFittingLogSolve::GetDescription() const
{
    std::string desc;

    desc = "Method templatefittinglogsolve:\n";

    desc.append("\tparam: templatefittinglogsolve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
    desc.append("\tparam: templatefittinglogsolve.overlapThreshold = <float value>\n");
    desc.append("\tparam: templatefittinglogsolve.extinction = {""yes"", ""no""}\n");
    desc.append("\tparam: templatefittinglogsolve.dustfit = {""yes"", ""no""}\n");
    desc.append("\tparam: templatefittinglogsolve.enablespclogrebin = {""yes"", ""no""}\n");
    desc.append("\tparam: templatefittinglogsolve.pdfcombination = {""marg"", ""bestchi2""}\n");
    desc.append("\tparam: templatefittinglogsolve.saveintermediateresults = {""yes"", ""no""}\n");

    return desc;

}


std::shared_ptr<CTemplateFittingSolveResult> CMethodTemplateFittingLogSolve::Compute(CDataStore& resultStore,
                                                                         const CSpectrum& spc,
                                                                         const CTemplateCatalog& tplCatalog,
                                                                         const TStringList& tplCategoryList,
                                                                         const TFloat64Range& lambdaRange,
                                                                         const TFloat64List& redshifts,
                                                                         Float64 overlapThreshold,
                                                                         std::vector<CMask> maskList,
                                                                         const string outputPdfRelDir,
                                                                         const Float64 redshiftSeparation,
                                                                         std::string spcComponent,
                                                                         std::string opt_interp,
                                                                         std::string opt_extinction,
                                                                         std::string opt_dustFit)
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "templatefittinglogsolve" );
    std::string scopeStr = "templatefitting";
    m_redshiftSeparation = redshiftSeparation;

    EType _type;
    if(spcComponent=="raw"){
       _type = nType_raw;
    }else if(spcComponent=="nocontinuum"){
       _type = nType_noContinuum;
       scopeStr = "templatefitting_nocontinuum";
    }else if(spcComponent=="continuum"){
        _type = nType_continuumOnly;
        scopeStr = "templatefitting_continuum";
    }else if(spcComponent=="all"){
        _type = nType_all;
    }

    resultStore.GetScopedParam( "extremacount", m_opt_maxCandidate, 5);
    resultStore.GetScopedParam( "pdfcombination", m_opt_pdfcombination, "marg");
    resultStore.GetScopedParam( "saveintermediateresults", m_opt_saveintermediateresults, "no");
    if(m_opt_saveintermediateresults=="yes")
    {
        m_opt_enableSaveIntermediateTemplateFittingResults = true;
    }else{
        m_opt_enableSaveIntermediateTemplateFittingResults = false;
    }
    resultStore.GetScopedParam( "enablespclogrebin", m_opt_spclogrebin, "yes");

    Log.LogInfo( "Method parameters:");
    Log.LogInfo( "    -overlapThreshold: %.3f", overlapThreshold);
    Log.LogInfo( "    -component: %s", spcComponent.c_str());
    Log.LogInfo( "    -IGM extinction: %s", opt_extinction.c_str());
    Log.LogInfo( "    -ISM dust-fit: %s", opt_dustFit.c_str());
    Log.LogInfo( "    -pdfcombination: %s", m_opt_pdfcombination.c_str());
    Log.LogInfo( "    -saveintermediateresults: %d", (int)m_opt_enableSaveIntermediateTemplateFittingResults);
    Log.LogInfo( "    -enable spectrum-log-rebin: %s", m_opt_spclogrebin.c_str());
    Log.LogInfo( "");

    if(m_opt_spclogrebin=="yes")
    {
        m_templateFittingOperator.enableSpcLogRebin(true);
    }else{
        m_templateFittingOperator.enableSpcLogRebin(false);
    }

    Log.LogInfo( "Iterating over %d tplCategories", tplCategoryList.size());
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

	Log.LogInfo( "   trying %s (%d templates)", category.c_str(), tplCatalog.GetTemplateCount( category ));
        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );

            Solve(resultStore, spc, tpl, lambdaRange, redshifts, overlapThreshold, maskList, _type, opt_interp, opt_extinction, opt_dustFit);

            storeResult = true;
        }
    }


    if( storeResult )
    {
        COperatorPdfz pdfz(m_opt_pdfcombination, m_redshiftSeparation, 0.0, m_opt_maxCandidate);

        std::shared_ptr<CPdfCandidateszResult> candidateResult = pdfz.Compute(BuildChisquareArray(resultStore, scopeStr));

        // save in resultstore pdf results
        std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
        resultStore.StoreGlobalResult( pdfPath.c_str(), pdfz.m_postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
        
        // save in resultstore candidates results
        {
            std::string name;
            if(resultStore.GetCurrentScopeName()=="templatefittinglogsolve")
                name = "candidatesresult";
            else  
                name = resultStore.GetCurrentScopeName() + "." +"candidatesresult";
            resultStore.StoreGlobalResult( name, candidateResult );
        }

        //for each candidate, get best model by reading from datastore and selecting best fit
        /////////////////////////////////////////////////////////////////////////////////////
        std::shared_ptr<const CExtremaResult> ExtremaResult = 
                    SaveExtremaResult( resultStore, scopeStr,
                                            candidateResult->m_ranked_candidates,
                                            spc,
                                            tplCatalog,
                                            tplCategoryList,
                                            lambdaRange,
                                            overlapThreshold,
                                            opt_interp,
                                            opt_extinction );

        // store extrema results
        StoreExtremaResults(resultStore, ExtremaResult);

        std::shared_ptr< CTemplateFittingSolveResult> TemplateFittingSolveResult = 
                        std::make_shared<CTemplateFittingSolveResult>(resultStore.GetCurrentScopeName(),
                                                                      ExtremaResult,
                                                                      m_opt_pdfcombination,
                                                                      pdfz.m_postmargZResult->valEvidenceLog);
        return TemplateFittingSolveResult;
    }

    return NULL;
}                                    


Bool CMethodTemplateFittingLogSolve::Solve(CDataStore& resultStore,
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
        //CRef<CTemplateFittingResult>  chisquareResult = (CTemplateFittingResult*)chiSquare.ExportChi2versusAZ( _spc, _tpl, lambdaRange, redshifts, overlapThreshold );
        auto  chisquareResult = std::dynamic_pointer_cast<CTemplateFittingResult>( m_templateFittingOperator.Compute( spc,
                                                                                                                       tpl,
                                                                                                                       lambdaRange,
                                                                                                                       redshifts,
                                                                                                                       overlapThreshold,
                                                                                                                       maskList,
                                                                                                                       opt_interp,
                                                                                                                       enable_extinction,
                                                                                                                       option_dustFitting ) );

        if( !chisquareResult )
        {
            //Log.LogError( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            resultStore.StoreScopedPerTemplateResult( tpl, scopeStr.c_str(), chisquareResult );

            //Save intermediate chisquare results
            if(m_opt_enableSaveIntermediateTemplateFittingResults && chisquareResult->ChiSquareIntermediate.size()>0 && chisquareResult->ChiSquareIntermediate.size()==chisquareResult->Redshifts.size())
            {
                Int32 nISM = chisquareResult->ChiSquareIntermediate[0].size();
                if(chisquareResult->ChiSquareIntermediate[0].size()>0)
                {
                    Int32 nIGM = chisquareResult->ChiSquareIntermediate[0][0].size();

                    for(Int32 kism=0; kism<nISM; kism++)
                    {
                        for(Int32 kigm=0; kigm<nIGM; kigm++)
                        {
                            std::shared_ptr<CTemplateFittingResult> result_chisquare_intermediate = std::shared_ptr<CTemplateFittingResult>( new CTemplateFittingResult() );
                            result_chisquare_intermediate->Init( chisquareResult->Redshifts.size(), 0, 0);
                            for(Int32 kz=0; kz<chisquareResult->Redshifts.size(); kz++)
                            {
                                result_chisquare_intermediate->Redshifts[kz] = chisquareResult->Redshifts[kz];
                                result_chisquare_intermediate->ChiSquare[kz] = chisquareResult->ChiSquareIntermediate[kz][kism][kigm];
                            }

                            std::string resname = (boost::format("%s_intermediate_ism%d_igm%d") % scopeStr.c_str() % kism % kigm).str();
                            resultStore.StoreScopedPerTemplateResult( tpl, resname.c_str(), result_chisquare_intermediate );
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

ChisquareArray CMethodTemplateFittingLogSolve::BuildChisquareArray(const CDataStore& store, const std::string & scopeStr) const
{
    ChisquareArray chisquarearray;

    Log.LogDetail("templatefittinglogsolve: build chisquare array");
    std::string scope = store.GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    templatefittinglogsolve: using results in scope: %s", scope.c_str());

    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    chisquarearray.cstLog = -1;

    Int32 retPdfz=-1;

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
            Log.LogInfo("templatefittinglogsolve: using cstLog = %f", chisquarearray.cstLog);
        }else if ( chisquarearray.cstLog != meritResult->CstLog)
        {
            Log.LogError("templatefittinglogsolve: Found different cstLog values in results... val-1=%f != val-2=%f", chisquarearray.cstLog, meritResult->CstLog);
            throw runtime_error("templatefittinglogsolve: Found different cstLog values in results");
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
                Log.LogError("templatefittinglogsolve: Found bad status result... for tpl=%s", (*it).first.c_str());
                throw runtime_error("templatefittinglogsolve: Found bad status result");
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
                Log.LogDetail("    templatefittinglogsolve: Pdfz combine - prepared merit  #%d for model : %s", chisquarearray.chisquares.size()-1, ((*it).first).c_str());
            }
        }
    }

    return chisquarearray;
}


std::shared_ptr<const CExtremaResult> 
CMethodTemplateFittingLogSolve::SaveExtremaResult(const CDataStore& store, const std::string & scopeStr,
                                               const TCandidateZbyRank & ranked_zCandidates,
                                               const CSpectrum& spc,
                                               const CTemplateCatalog& tplCatalog,
                                               const TStringList& tplCategoryList,
                                               const TFloat64Range& lambdaRange,
                                               Float64 overlapThreshold,
                                               std::string opt_interp,
                                               std::string opt_extinction
                                               )
{

    Log.LogDetail("CMethodTemplateFittingLogSolve::SaveExtremaResult: building chisquare array");
    std::string scope = store.GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    templatefittinglogsolve: using results in scope: %s", scope.c_str());

    TOperatorResultMap results = store.GetPerTemplateResult(scope.c_str());

    //Bool foundRedshiftAtLeastOnce = false;

    Int32 extremumCount = ranked_zCandidates.size();
    
    auto firstResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( (*results.begin()).second );
    const TFloat64List & redshifts = firstResult->Redshifts;

    //check all results  status
    for ( auto & r : results)
    {
        auto TplFitResult = std::dynamic_pointer_cast<const CTemplateFittingResult>( r.second );
        if(TplFitResult->ChiSquare.size() != redshifts.size()){
            Log.LogError("CMethodTemplateFittingLogSolve::SaveExtremaResult, templatefitting results (for tpl=%s) has wrong size", r.first.c_str());
            throw runtime_error("CMethodTemplateFittingLogSolve::SaveExtremaResult, one templatefitting results has wrong size");
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
            Log.LogError("CMethodTemplateFittingLogSolve::SaveExtremaResult: Found bad status result... for tpl=%s", r.first.c_str());
            throw runtime_error("CMethodTemplateFittingLogSolve::SaveExtremaResult: Found bad status result");
        }
        if(foundBadRedshift)
        {
            Log.LogError("CMethodTemplateFittingLogSolve::SaveExtremaResult: redshift vector is not the same for tpl=%s", r.first.c_str());
            throw runtime_error("CMethodTemplateFittingLogSolve::SaveExtremaResult: Found different redshift vector");
        }

    }

    std::shared_ptr<CExtremaResult> ExtremaResult = make_shared<CExtremaResult>(extremumCount);

    ExtremaResult->Candidates = ranked_zCandidates;

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
        ExtremaResult->FittedTplMerit[i] = ChiSquare;
        ExtremaResult->FittedTplName[i] = tplName;
        ExtremaResult->FittedTplMeiksinIdx[i] = TplFitResult->FitMeiksinIdx[idx];
        ExtremaResult->FittedTplDustCoeff[i] = TplFitResult->FitDustCoeff[idx];
        ExtremaResult->FittedTplAmplitude[i] = TplFitResult->FitAmplitude[idx];
        ExtremaResult->FittedTplAmplitudeError[i] = TplFitResult->FitAmplitudeError[idx];
        ExtremaResult->FittedTplDtm[i] = TplFitResult->FitDtM[idx];
        ExtremaResult->FittedTplMtm[i] = TplFitResult->FitMtM[idx];
        ExtremaResult->FittedTplLogPrior[i] = TplFitResult->LogPrior[idx];

        Float64 FitSNR = NAN;
        if (TplFitResult->FitMtM[idx] != 0.)
            FitSNR = abs(TplFitResult->FitDtM[idx])/sqrt(TplFitResult->FitMtM[idx]); // = |amplitude|/amplitudeError
        ExtremaResult->FittedTplSNR[i] = FitSNR;

        const CTemplate& tpl = tplCatalog.GetTemplateByName(tplCategoryList, tplName);
        std::shared_ptr<CModelSpectrumResult> spcmodelPtr; 

        COperatorTemplateFitting templateFittingOperator;
        templateFittingOperator.ComputeSpectrumModel(spc, tpl, 
                                                        z,
                                                        TplFitResult->FitDustCoeff[idx],
                                                        TplFitResult->FitMeiksinIdx[idx],
                                                        TplFitResult->FitAmplitude[idx],
                                                        opt_interp, opt_extinction, lambdaRange, 
                                                        overlapThreshold, spcmodelPtr);
        ExtremaResult->m_savedModelSpectrumResults[i] = std::move(spcmodelPtr);

        ExtremaResult->m_savedModelContinuumFittingResults[i] = 
                        std::make_shared<CModelContinuumFittingResult>( z,
                                                                        tplName,
                                                                        ChiSquare,
                                                                        TplFitResult->FitAmplitude[idx],
                                                                        TplFitResult->FitAmplitudeError[idx],
                                                                        TplFitResult->FitDustCoeff[idx],
                                                                        TplFitResult->FitMeiksinIdx[idx],
                                                                        FitSNR );          
    }

    return ExtremaResult;
}


void CMethodTemplateFittingLogSolve::StoreExtremaResults( CDataStore &dataStore, 
                                                       std::shared_ptr<const CExtremaResult> & ExtremaResult) const
{
    Int32 nResults = ExtremaResult->size();

    auto & m_savedModelContinuumFittingResults = ExtremaResult->m_savedModelContinuumFittingResults;
    auto & m_savedModelSpectrumResults = ExtremaResult->m_savedModelSpectrumResults;

    if( m_savedModelContinuumFittingResults.size()!= m_savedModelSpectrumResults.size()){
        Log.LogError(" CMethodTemplateFittingLogSolve::SaveSpectrumResults: spectrumModel size doesnt not correspond to modelParam size.");
        throw runtime_error(" CMethodTemplateFittingLogSolve::SaveSpectrumResults: spectrumModel size doesnt not correspond to modelParam size. Aborting!");
    }

    Log.LogDetail("  Methode-COperatorTemplatefittingLog: now saving spectrum/model templatefitting results n=%d", m_savedModelSpectrumResults.size());
    for(Int32 i=0; i<m_savedModelSpectrumResults.size(); i++)
    {
        std::string fname_spc = (boost::format("templatefitting_spc_extrema_%1%") % i).str();
        dataStore.StoreScopedGlobalResult( fname_spc.c_str(), m_savedModelSpectrumResults[i] );

        fname_spc = (boost::format("templatefitting_fitcontinuum_extrema_%1%") % i).str();
        dataStore.StoreScopedGlobalResult( fname_spc.c_str(),  m_savedModelContinuumFittingResults[i] );
    }
    return;
}

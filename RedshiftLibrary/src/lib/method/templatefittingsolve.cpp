#include <RedshiftLibrary/method/templatefittingsolve.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/statistics/zprior.h>
#include <RedshiftLibrary/statistics/deltaz.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <cfloat>

using namespace NSEpic;
using namespace std;


CMethodTemplateFittingSolve::CMethodTemplateFittingSolve()
{
}

CMethodTemplateFittingSolve::~CMethodTemplateFittingSolve()
{
}

const std::string CMethodTemplateFittingSolve::GetDescription()
{
    std::string desc;

    desc = "Method templatefittingsolve:\n";

    desc.append("\tparam: templatefittingsolve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
    desc.append("\tparam: templatefittingsolve.overlapThreshold = <float value>\n");
    desc.append("\tparam: templatefittingsolve.interpolation = {""precomputedfinegrid"", ""lin""}\n");
    desc.append("\tparam: templatefittingsolve.extinction = {""yes"", ""no""}\n");
    desc.append("\tparam: templatefittingsolve.dustfit = {""yes"", ""no""}\n");
    desc.append("\tparam: templatefittingsolve.pdfcombination = {""marg"", ""bestchi2""}\n");
    desc.append("\tparam: templatefittingsolve.saveintermediateresults = {""yes"", ""no""}\n");

    return desc;

}


std::shared_ptr<CTemplateFittingSolveResult> CMethodTemplateFittingSolve::Compute(CDataStore& resultStore,
                                                                       const CSpectrum& spc,
                                                                       const CTemplateCatalog& tplCatalog,
                                                                       const TStringList& tplCategoryList,
                                                                       const TFloat64Range& lambdaRange,
                                                                       const TFloat64List& redshifts,
                                                                       Float64 overlapThreshold,
                                                                       std::vector<CMask> maskList,
                                                                       const std::string outputPdfRelDir,
                                                                       const Float64 radius,
                                                                       std::string spcComponent,
                                                                       std::string opt_interp,
                                                                       std::string opt_extinction,
                                                                       std::string opt_dustFit)
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "templatefittingsolve" );
    std::string scopeStr = "templatefitting";
    m_radius = radius;

    CTemplateFittingSolveResult::EType _type;
    if(spcComponent=="raw"){
       _type = CTemplateFittingSolveResult::nType_raw;
    }else if(spcComponent=="nocontinuum"){
       _type = CTemplateFittingSolveResult::nType_noContinuum;
       scopeStr = "templatefitting_nocontinuum";
    }else if(spcComponent=="continuum"){
        _type = CTemplateFittingSolveResult::nType_continuumOnly;
        scopeStr = "templatefitting_continuum";
    }else if(spcComponent=="all"){
        _type = CTemplateFittingSolveResult::nType_all;
    }

    resultStore.GetScopedParam( "pdfcombination", m_opt_pdfcombination, "marg");
    resultStore.GetScopedParam( "saveintermediateresults", m_opt_saveintermediateresults, "no");
    if(m_opt_saveintermediateresults=="yes")
    {
        m_opt_enableSaveIntermediateTemplateFittingResults = true;
    }else{
        m_opt_enableSaveIntermediateTemplateFittingResults = false;
    }

    Log.LogInfo( "Method parameters:");
    Log.LogInfo( "    -overlapThreshold: %.3f", overlapThreshold);
    Log.LogInfo( "    -component: %s", spcComponent.c_str());
    Log.LogInfo( "    -interp: %s", opt_interp.c_str());
    Log.LogInfo( "    -IGM extinction: %s", opt_extinction.c_str());
    Log.LogInfo( "    -ISM dust-fit: %s", opt_dustFit.c_str());
    Log.LogInfo( "    -pdfcombination: %s", m_opt_pdfcombination.c_str());
    Log.LogInfo( "    -saveintermediateresults: %d", (int)m_opt_enableSaveIntermediateTemplateFittingResults);
    Log.LogInfo( "");

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
        std::shared_ptr< CTemplateFittingSolveResult> TemplateFittingSolveResult = std::shared_ptr< CTemplateFittingSolveResult>( new CTemplateFittingSolveResult(_type, resultStore.GetCurrentScopeName()) );

        std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
        Int32 retCombinePdf = CombinePDF(resultStore, scopeStr, m_opt_pdfcombination, postmargZResult);

        if(retCombinePdf==0)
        {
            //check pdf sum=1
            CPdfz pdfz;
            Float64 sumRect = pdfz.getSumRect(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Float64 sumTrapez = pdfz.getSumTrapez(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Log.LogDetail("    templatefittingsolve: Pdfz normalization - sum rect. = %e", sumRect);
            Log.LogDetail("    templatefittingsolve: Pdfz normalization - sum trapz. = %e", sumTrapez);
            Bool pdfSumCheck = abs(sumRect-1.0)<1e-1 || abs(sumTrapez-1.0)<1e-1;
            if(!pdfSumCheck){
                Log.LogError("    templatefittingsolve: Pdfz normalization failed (rectsum = %f, trapzesum = %f)", sumRect, sumTrapez);
            }


            {
                std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
                resultStore.StoreGlobalResult( pdfPath.c_str(), postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
            }
        }
        Int32 n_cand = 5; //this is hardcoded for now for this method
        std::vector<Float64> zcandidates_unordered_list;
        Bool retzc = TemplateFittingSolveResult->GetRedshiftCandidates( resultStore, zcandidates_unordered_list, n_cand, outputPdfRelDir.c_str());
        if(retzc)
        {
                Log.LogInfo( "Found %d z-candidates", zcandidates_unordered_list.size() );
        }else{
                Log.LogError( "Failed to get z candidates from these results");
        }

        Bool b = ExtractCandidateResults(resultStore, zcandidates_unordered_list, outputPdfRelDir.c_str());
        //for each candidate, get best model by reading from resultstore and selecting best fit
        for(Int32 i = 0; i<zcandidates_unordered_list.size(); i++){
            Int32 ret = TemplateFittingSolveResult->GetBestModel(resultStore, zcandidates_unordered_list[i]); //, tplName, MeiksinIdx, DustCoeff, Amplitude);

            if(ret==-1){
                Log.LogError("  TemplateFittingSolve: Couldn't find best model for candidate %f", zcandidates_unordered_list[i]);
                continue;
            }
            //now that we have best tplName, we have access to meiksin index & EBMV to create the model spectrum
            try {
                const CTemplate& tpl = tplCatalog.GetTemplateByName(tplCategoryList, TemplateFittingSolveResult->GetTemplateName());
                CModelSpectrumResult spcmodel;
                m_templateFittingOperator.ComputeSpectrumModel(spc, tpl, 
                                                 zcandidates_unordered_list[i],
                                                 TemplateFittingSolveResult->GetDustCoeff(), 
                                                 TemplateFittingSolveResult->GetMeiksinIdx(), 
                                                 TemplateFittingSolveResult->GetAmplitude(),
                                                 opt_interp, opt_extinction, lambdaRange, 
                                                 overlapThreshold, spcmodel);
                m_savedModelSpectrumResults.push_back(std::make_shared<CModelSpectrumResult>(spcmodel));
                m_savedModelContinuumFittingResults.push_back(std::make_shared<CModelContinuumFittingResult>(
                                                                                zcandidates_unordered_list[i],
                                                                                tpl.GetName(), 
                                                                                TemplateFittingSolveResult->GetMerit(),
                                                                                TemplateFittingSolveResult->GetAmplitude(),
                                                                                TemplateFittingSolveResult->GetAmplitudeError(),
                                                                                TemplateFittingSolveResult->GetDustCoeff(),
                                                                                TemplateFittingSolveResult->GetMeiksinIdx(),
                                                                                TemplateFittingSolveResult->GetFittingSNR()));          
            }catch(const std::runtime_error& e){
                Log.LogError("  TemplateFittingSolve: Couldn't find template by tplName: %s for candidate %f", TemplateFittingSolveResult->GetTemplateName().c_str(), zcandidates_unordered_list[i]);
                continue;
            }   
                    
        }
        SaveSpectrumResults(resultStore);
        return TemplateFittingSolveResult;
    }

    return NULL;
}

Bool CMethodTemplateFittingSolve::Solve(CDataStore& resultStore,
                                   const CSpectrum& spc,
                                   const CTemplate& tpl,
                                   const TFloat64Range& lambdaRange,
                                   const TFloat64List& redshifts,
                                   Float64 overlapThreshold,
                                   std::vector<CMask> maskList,
                                   CTemplateFittingSolveResult::EType spctype,
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
    if(spctype == CTemplateFittingSolveResult::nType_all){
        _ntype = 3;
    }

    const CSpectrum::EType save_spcType = spc.GetType();
    const CSpectrum::EType save_tplType = tpl.GetType();

    for( Int32 i=0; i<_ntype; i++){
        if(spctype == CTemplateFittingSolveResult::nType_all){
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
        auto  templateFittingResult = std::dynamic_pointer_cast<CTemplateFittingResult>( m_templateFittingOperator.Compute( spc,
                                                                                                           tpl,
                                                                                                           lambdaRange,
                                                                                                           redshifts,
                                                                                                           overlapThreshold,
                                                                                                           maskList,
                                                                                                           opt_interp,
                                                                                                           enable_extinction,
                                                                                                           option_dustFitting ) );

        templateFittingResult->CallFindExtrema(m_radius);
        
        if( !templateFittingResult )
        {
            //Log.LogError( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            resultStore.StoreScopedPerTemplateResult( tpl, scopeStr.c_str(), templateFittingResult );

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

Int32 CMethodTemplateFittingSolve::CombinePDF(CDataStore& store, std::string scopeStr, std::string opt_combine, std::shared_ptr<CPdfMargZLogResult> postmargZResult)
{
    Log.LogInfo("templatefittingsolve: Pdfz computation");
    std::string scope = store.GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    templatefittingsolve: using results in scope: %s", scope.c_str());

    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    CPdfz pdfz;
    Float64 cstLog = -1;

    Int32 retPdfz=-1;
    std::vector<TFloat64List> priors;
    std::vector<TFloat64List> chiSquares;
    std::vector<Float64> redshifts;
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
        if(cstLog==-1)
        {
            cstLog = meritResult->CstLog;
            Log.LogInfo("templatefittingsolve: using cstLog = %f", cstLog);
        }else if ( cstLog != meritResult->CstLog)
        {
            Log.LogError("templatefittingsolve: Found different cstLog values in results... val-1=%f != val-2=%f", cstLog, meritResult->CstLog);
        }
        if(redshifts.size()==0)
        {
            redshifts = meritResult->Redshifts;
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
            }
        }

        CZPrior zpriorhelper;
        for(Int32 kism=0; kism<nISM; kism++)
        {
            for(Int32 kigm=0; kigm<nIGM; kigm++)
            {
                priors.push_back(zpriorhelper.GetConstantLogZPrior(meritResult->Redshifts.size()));

                //correct chi2 for ampl. marg. if necessary: todo add switch, currently deactivated
                chiSquares.emplace_back(meritResult->ChiSquareIntermediate.size(), DBL_MAX);
                TFloat64List & logLikelihoodCorrected= chiSquares.back();
                for ( UInt32 kz=0; kz<meritResult->Redshifts.size(); kz++)
                {
                    logLikelihoodCorrected[kz] = meritResult->ChiSquareIntermediate[kz][kism][kigm];// + resultXXX->ScaleMargCorrectionTplshapes[][]?;
                }
                Log.LogDetail("    templatefittingsolve: Pdfz combine - prepared merit  #%d for model : %s", chiSquares.size()-1, ((*it).first).c_str());
            }
        }
    }

    if(chiSquares.size()>0)
    {
        if(opt_combine=="marg")
        {
            Log.LogInfo("    templatefittingsolve: Pdfz combination - Marginalization");
            retPdfz = pdfz.Marginalize( redshifts, chiSquares, priors, cstLog, postmargZResult);
        }else if(opt_combine=="bestchi2")
        {
            Log.LogInfo("    templatefittingsolve: Pdfz combination - BestChi2");
            retPdfz = pdfz.BestChi2( redshifts, chiSquares, priors, cstLog, postmargZResult);
        }else{
            Log.LogError("    templatefittingsolve: Unable to parse pdf combination method option");
        }
    }else
    {
        Log.LogError("    templatefittingsolve: Unable to find any chisquares prepared for combination. chiSquares.size()=%d", chiSquares.size());
    }


    if(retPdfz!=0)
    {
        Log.LogError("    templatefittingsolve: Pdfz computation failed");
    }

    return retPdfz;
}

Bool CMethodTemplateFittingSolve::ExtractCandidateResults(CDataStore &store, std::vector<Float64> zcandidates_unordered_list, std::string outputPdfRelDir)
{
        Log.LogInfo( "Computing candidates Probabilities" );
        std::shared_ptr<CPdfCandidateszResult> zcand = std::shared_ptr<CPdfCandidateszResult>(new CPdfCandidateszResult());

        std::string scope_res = outputPdfRelDir + "/logposterior.logMargP_Z_data";
        auto results =  store.GetGlobalResult( scope_res.c_str() );
        auto logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results.lock() );

        if(!logzpdf1d)
        {
            Log.LogError( "Extract Proba. for z candidates: no results retrieved from scope: %s", scope_res.c_str());
            throw std::runtime_error("Extract Proba. for z candidates: no results retrieved from scope");
        }

        //Compute Deltaz should happen after marginalization
        // use it for computing the integrated PDF
        TFloat64List deltaz;
        CDeltaz deltaz_obj;
        for(Int32 i =0; i<zcandidates_unordered_list.size(); i++){
            Float64 z = zcandidates_unordered_list[i];
            deltaz.push_back(deltaz_obj.GetDeltaz(logzpdf1d->Redshifts, logzpdf1d->valProbaLog, z));
        }

        Log.LogInfo( "  Integrating %d candidates proba.", zcandidates_unordered_list.size() );
        zcand->Compute(zcandidates_unordered_list, logzpdf1d->Redshifts, logzpdf1d->valProbaLog, deltaz);
         std::string name;
        if(store.GetCurrentScopeName()=="templatefittingsolve")
            name = "candidatesresult";
        else  
            name = store.GetCurrentScopeName() + "." +"candidatesresult";
        store.StoreGlobalResult( name, zcand ); 

    return true;
}

void CMethodTemplateFittingSolve::SaveSpectrumResults(CDataStore &dataStore)
{
    if( m_savedModelContinuumFittingResults.size()!= m_savedModelSpectrumResults.size()){
        Log.LogError(" CMethodTemplateFittingSolve::SaveSpectrumResults: spectrumModel size doesnt not correspond to modelParam size.");
        throw runtime_error(" CMethodTemplateFittingSolve::SaveSpectrumResults: spectrumModel size doesnt not correspond to modelParam size. Aborting!");
    }

    Log.LogDetail("  Methode-COperatorTemplatefitting: now saving spectrum/model templatefitting results n=%d", m_savedModelSpectrumResults.size());
    for(Int32 i=0; i<m_savedModelSpectrumResults.size(); i++)
    {
        std::string fname_spc = (boost::format("templatefitting_spc_extrema_%1%") % i).str();
        dataStore.StoreScopedGlobalResult( fname_spc.c_str(), m_savedModelSpectrumResults[i] );

        fname_spc = (boost::format("templatefitting_fitcontinuum_extrema_%1%") % i).str();
        dataStore.StoreScopedGlobalResult( fname_spc.c_str(),  m_savedModelContinuumFittingResults[i] );
    }
    return;
}

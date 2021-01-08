#include <RedshiftLibrary/method/templatefittinglogsolve.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/statistics/deltaz.h>
#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/operator/pdfLogresult.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <RedshiftLibrary/common/quicksort.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <cfloat>

using namespace NSEpic;
using namespace std;


CMethodTemplateFittingLogSolve::CMethodTemplateFittingLogSolve( std::string calibrationPath )
{
    m_templateFittingOperator = new COperatorTemplateFittingLog( calibrationPath );
}

CMethodTemplateFittingLogSolve::~CMethodTemplateFittingLogSolve()
{
    delete m_templateFittingOperator;
}

const std::string CMethodTemplateFittingLogSolve::GetDescription()
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
                                                                         const Float64 radius,
                                                                         std::string spcComponent,
                                                                         std::string opt_interp,
                                                                         std::string opt_extinction,
                                                                         std::string opt_dustFit)
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "templatefittinglogsolve" );
    std::string scopeStr = "chisquare";
    m_radius = radius;

    CTemplateFittingSolveResult::EType _type;
    if(spcComponent=="raw"){
       _type = CTemplateFittingSolveResult::nType_raw;
    }else if(spcComponent=="nocontinuum"){
       _type = CTemplateFittingSolveResult::nType_noContinuum;
       scopeStr = "chisquare_nocontinuum";
    }else if(spcComponent=="continuum"){
        _type = CTemplateFittingSolveResult::nType_continuumOnly;
        scopeStr = "chisquare_continuum";
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
        m_templateFittingOperator->enableSpcLogRebin(true);
    }else{
        m_templateFittingOperator->enableSpcLogRebin(false);
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
        std::shared_ptr< CTemplateFittingSolveResult> TemplateFittingSolveResult = std::shared_ptr< CTemplateFittingSolveResult>( new CTemplateFittingSolveResult(_type, "templatefittinglogsolve") );

        std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
        Int32 retCombinePdf = CombinePDF(resultStore, scopeStr, m_opt_pdfcombination, postmargZResult);

        if(retCombinePdf==0)
        {
            //check pdf sum=1
            CPdfz pdfz;
            Float64 sumRect = pdfz.getSumRect(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Float64 sumTrapez = pdfz.getSumTrapez(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Log.LogDetail("    templatefittinglogsolve: Pdfz normalization - sum rect. = %e", sumRect);
            Log.LogDetail("    templatefittinglogsolve: Pdfz normalization - sum trapz. = %e", sumTrapez);
            Bool pdfSumCheck = abs(sumRect-1.0)<1e-1 || abs(sumTrapez-1.0)<1e-1;
            if(!pdfSumCheck){
                Log.LogError("    templatefittinglogsolve: Pdfz normalization failed (rectsum = %f, trapzesum = %f)", sumRect, sumTrapez);
            }


            {
                std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
                resultStore.StoreGlobalResult( pdfPath.c_str(), postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
            }
        }

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
                                     CTemplateFittingSolveResult::EType spctype,
                                     std::string opt_interp,
                                     std::string opt_extinction,
                                     std::string opt_dustFitting)
{
    std::string scopeStr = "chisquare";
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
            scopeStr = "chisquare_continuum";

        }else if(_spctype == CSpectrum::nType_raw){
            // use full spectrum
            scopeStr = "chisquare";

        }else if(_spctype == CSpectrum::nType_noContinuum){
            // use spectrum without continuum
            scopeStr = "chisquare_nocontinuum";
            //
            option_dustFitting = -1;
        }

        // Compute merit function
        //CRef<CTemplateFittingResult>  chisquareResult = (CTemplateFittingResult*)chiSquare.ExportChi2versusAZ( _spc, _tpl, lambdaRange, redshifts, overlapThreshold );
        auto  chisquareResult = std::dynamic_pointer_cast<CTemplateFittingResult>( m_templateFittingOperator->Compute( spc,
                                                                                                                       tpl,
                                                                                                                       lambdaRange,
                                                                                                                       redshifts,
                                                                                                                       overlapThreshold,
                                                                                                                       maskList,
                                                                                                                       opt_interp,
                                                                                                                       enable_extinction,
                                                                                                                       option_dustFitting ) );
        chisquareResult->CallFindExtrema(m_radius);

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

Int32 CMethodTemplateFittingLogSolve::CombinePDF(CDataStore& store, std::string scopeStr, std::string opt_combine, std::shared_ptr<CPdfMargZLogResult> postmargZResult)
{
    Log.LogInfo("templatefittinglogsolve: Pdfz computation");
    std::string scope = store.GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    templatefittinglogsolve: using results in scope: %s", scope.c_str());

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
            Log.LogInfo("templatefittinglogsolve: using cstLog = %f", cstLog);
        }else if ( cstLog != meritResult->CstLog)
        {
            Log.LogError("templatefittinglogsolve: Found different cstLog values in results... val-1=%f != val-2=%f", cstLog, meritResult->CstLog);
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
                Log.LogError("templatefittinglogsolve: Found bad status result... for tpl=%s", (*it).first.c_str());
            }
        }

        for(Int32 kism=0; kism<nISM; kism++)
        {
            for(Int32 kigm=0; kigm<nIGM; kigm++)
            {
                priors.push_back(pdfz.GetConstantLogZPrior(meritResult->Redshifts.size()));

                //correct chi2 for ampl. marg. if necessary: todo add switch, currently deactivated
                chiSquares.emplace_back(meritResult->ChiSquareIntermediate.size(), DBL_MAX);
                TFloat64List & logLikelihoodCorrected = chiSquares.back();
                for ( UInt32 kz=0; kz<meritResult->Redshifts.size(); kz++)
                {
                    logLikelihoodCorrected[kz] = meritResult->ChiSquareIntermediate[kz][kism][kigm];// + resultXXX->ScaleMargCorrectionTplshapes[][]?;
                }
                Log.LogDetail("    templatefittinglogsolve: Pdfz combine - prepared merit  #%d for model : %s, ism=%d, igm=%d", chiSquares.size()-1, ((*it).first).c_str(), kism, kigm);
            }
        }
    }

    if(chiSquares.size()>0)
    {
        if(opt_combine=="marg")
        {
            Log.LogInfo("    templatefittinglogsolve: Pdfz combination - Marginalization");
            retPdfz = pdfz.Marginalize( redshifts, chiSquares, priors, cstLog, postmargZResult);
        }else if(opt_combine=="bestchi2")
        {
            Log.LogInfo("    templatefittinglogsolve: Pdfz combination - BestChi2");
            retPdfz = pdfz.BestChi2( redshifts, chiSquares, priors, cstLog, postmargZResult);
        }else{
            Log.LogError("    templatefittinglogsolve: Unable to parse pdf combination method option");
        }
    }else
    {
        Log.LogError("    templatefittinglogsolve: Unable to find any chisquares prepared for combination. chiSquares.size()=%d", chiSquares.size());
    }


    if(retPdfz!=0)
    {
        Log.LogError("    templatefittinglogsolve: Pdfz computation failed");
    }

    return retPdfz;
}

Bool CMethodTemplateFittingLogSolve::ExtractCandidateResults(CDataStore& store, std::vector<Float64> zcandidates_unordered_list)
{
        Log.LogInfo( "Computing candidates Probabilities" );
        std::shared_ptr<CPdfCandidateszResult> zcand = std::shared_ptr<CPdfCandidateszResult>(new CPdfCandidateszResult());

        std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
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
        
        store.StoreScopedGlobalResult( "candidatesresult", zcand ); 

    return true;
}

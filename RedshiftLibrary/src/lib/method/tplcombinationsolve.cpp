#include <RedshiftLibrary/method/tplcombinationsolve.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <RedshiftLibrary/statistics/deltaz.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <float.h>
using namespace NSEpic;
using namespace std;

CMethodTplcombinationSolve::CMethodTplcombinationSolve( std::string calibrationPath )
{
    m_tplcombinationOperator = new COperatorTplcombination( calibrationPath );
}

CMethodTplcombinationSolve::~CMethodTplcombinationSolve()
{
    delete m_tplcombinationOperator;
}


const std::string CMethodTplcombinationSolve::GetDescription()
{
    std::string desc;

    desc = "Method tplcombinationsolve:\n";

    desc.append("\tparam: tplcombinationsolve.spectrum.component = {""raw"", ""nocontinuum"", ""continuum"", ""all""}\n");
    desc.append("\tparam: tplcombinationsolve.overlapThreshold = <float value>\n");
    desc.append("\tparam: tplcombinationsolve.interpolation = {""precomputedfinegrid"", ""lin""}\n");
    //desc.append("\tparam: tplcombinationsolve.extinction = {""yes"", ""no""}\n");
    //desc.append("\tparam: tplcombinationsolve.dustfit = {""yes"", ""no""}\n");
    //desc.append("\tparam: tplcombinationsolve.pdfcombination = {""marg"", ""bestchi2""}\n");
    desc.append("\tparam: tplcombinationsolve.saveintermediateresults = {""yes"", ""no""}\n");


    return desc;

}


std::shared_ptr<CChisquareSolveResult> CMethodTplcombinationSolve::Compute(CDataStore& resultStore,
                                                                              const CSpectrum& spc,
                                                                              const CSpectrum& spcWithoutCont,
                                                                              const CTemplateCatalog& tplCatalog,
                                                                              const TStringList& tplCategoryList,
                                                                              const TFloat64Range& lambdaRange,
                                                                              const TFloat64List& redshifts,
                                                                              Float64 overlapThreshold,
                                                                              std::vector<CMask> maskList,
                                                                              const std::string outputPdfRelDir,
                                                                              std::string spcComponent,
                                                                              std::string opt_interp,
                                                                              std::string opt_extinction,
                                                                              std::string opt_dustFit)
{
    Bool storeResult = false;

//    std::string _name = "Tplcombination";
    CDataStore::CAutoScope resultScope( resultStore, "tplcombinationsolve" );
//    std::string _scope = "tplcombination";
    std::string scopeStr = "chisquare";

    Int32 _type;
    if(spcComponent=="raw"){
       _type = CChisquareSolveResult::nType_raw;
    }else if(spcComponent=="nocontinuum"){
       _type = CChisquareSolveResult::nType_noContinuum;
       scopeStr = "chisquare_nocontinuum";
    }else if(spcComponent=="continuum"){
        _type = CChisquareSolveResult::nType_continuumOnly;
        scopeStr = "chisquare_continuum";
    }else if(spcComponent=="all"){
        _type = CChisquareSolveResult::nType_all;
    }


    resultStore.GetScopedParam( "pdfcombination", m_opt_pdfcombination, "marg");
    resultStore.GetScopedParam( "saveintermediateresults", m_opt_saveintermediateresults, "no");
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
    Log.LogInfo( "    -component: %s", spcComponent.c_str());
    Log.LogInfo( "    -IGM extinction: %s", opt_extinction.c_str());
    Log.LogInfo( "    -ISM dust-fit: %s", opt_dustFit.c_str());
    //Log.LogInfo( "    -pdfcombination: %s", m_opt_pdfcombination.c_str());
    Log.LogInfo( "    -saveintermediateresults: %d", (int)m_opt_enableSaveIntermediateChisquareResults);
    Log.LogInfo( "");

//    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
//    {
//        std::string category = tplCategoryList[i];

//        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
//        {
//            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
//            const CTemplate& tplWithoutCont = tplCatalog.GetTemplateWithoutContinuum( category, j );

    Solve( resultStore, spc, spcWithoutCont, tplCatalog, tplCategoryList, lambdaRange, redshifts, overlapThreshold, maskList, _type, opt_interp, opt_extinction, opt_dustFit);

    storeResult = true;
//        }
//    }


    if( storeResult )
    {
        std::shared_ptr< CChisquareSolveResult>  solveResult =
                std::shared_ptr< CChisquareSolveResult>( new CChisquareSolveResult(_type, "tplcombinationsolve") );

        std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
        Int32 retCombinePdf = CombinePDF(resultStore, scopeStr, m_opt_pdfcombination, postmargZResult);

        if(retCombinePdf==0)
        {
            //check pdf sum=1
            CPdfz pdfz;
            Float64 sumRect = pdfz.getSumRect(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Float64 sumTrapez = pdfz.getSumTrapez(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Log.LogDetail("    tplcombinationsolve: Pdfz normalization - sum rect. = %e", sumRect);
            Log.LogDetail("    tplcombinationsolve: Pdfz normalization - sum trapz. = %e", sumTrapez);
            Bool pdfSumCheck = abs(sumRect-1.0)<1e-1 || abs(sumTrapez-1.0)<1e-1;
            if(!pdfSumCheck){
                Log.LogError("    tplcombinationsolve: Pdfz normalization failed (rectsum = %f, trapzesum = %f)", sumRect, sumTrapez);
            }


            {
                std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
                resultStore.StoreGlobalResult( pdfPath.c_str(), postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
            }
        }

        return solveResult;
    }

    return NULL;
}

Bool CMethodTplcombinationSolve::Solve(CDataStore& resultStore,
                                    const CSpectrum& spc,
                                    const CSpectrum& spcWithoutCont,
                                    const CTemplateCatalog &tplCatalog,
                                    const TStringList &tplCategoryList,
                                    const TFloat64Range& lambdaRange,
                                    const TFloat64List& redshifts,
                                    Float64 overlapThreshold,
                                    std::vector<CMask> maskList,
                                    Int32 spctype,
                                    std::string opt_interp,
                                    std::string opt_extinction,
                                    std::string opt_dustFitting )
{
    CSpectrum _spc;
    std::string scopeStr = "tplcombination";
    Int32 _ntype = 1;
    Int32 _spctype = spctype;
    Int32 _spctypetab[3] = {CChisquareSolveResult::nType_raw, CChisquareSolveResult::nType_noContinuum, CChisquareSolveResult::nType_continuumOnly};


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
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            CTemplate _tpl=tpl;
            tplList.push_back(_tpl);
        }
    }

    //case: nType_all
    if(spctype == CChisquareSolveResult::nType_all){
        _ntype = 3;
    }
    for( Int32 i=0; i<_ntype; i++){
        if(spctype == CChisquareSolveResult::nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = spctype;
        }

        if(_spctype == CChisquareSolveResult::nType_continuumOnly){
            // use continuum only
            _spc = spc;
            CSpectrumFluxAxis spcfluxAxis = _spc.GetFluxAxis();
            spcfluxAxis.Subtract(spcWithoutCont.GetFluxAxis());
            CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
            sfluxAxisPtr = spcfluxAxis;
            scopeStr = "chisquare_continuum";
        }else if(_spctype == CChisquareSolveResult::nType_raw){
            // use full spectrum
            _spc = spc;
            scopeStr = "chisquare";

        }else if(_spctype == CChisquareSolveResult::nType_noContinuum){
            // use spectrum without continuum
            _spc = spc;
            CSpectrumFluxAxis spcfluxAxis = spcWithoutCont.GetFluxAxis();
            CSpectrumFluxAxis& sfluxAxisPtr = _spc.GetFluxAxis();
            sfluxAxisPtr = spcfluxAxis;
            scopeStr = "chisquare_nocontinuum";
            //
            enable_dustFitting = 0;
        }

        // Compute merit function
        auto  result = std::dynamic_pointer_cast<CChisquareResult>( m_tplcombinationOperator->Compute( _spc, tplList, lambdaRange, redshifts, overlapThreshold, maskList, opt_interp, enable_extinction, enable_dustFitting ) );

        if( !result )
        {
            //Log.LogInfo( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            Log.LogDetail("tplcombinationsolve: Save tplcombination results");
            resultStore.StoreScopedGlobalResult(scopeStr.c_str(), result );
            // Store spectrum results
            Log.LogDetail("tplcombinationsolve: Save spectrum/model results");
            m_tplcombinationOperator->SaveSpectrumResults(resultStore);
        }
    }

    return true;
}


Int32 CMethodTplcombinationSolve::CombinePDF(CDataStore &store, std::string scopeStr, std::string opt_combine, std::shared_ptr<CPdfMargZLogResult> postmargZResult)
{
    Log.LogInfo("tplcombinationsolve: Pdfz computation");
    std::string scope = store.GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDebug("tplcombinationsolve: combining pdf using results at scope: %s", scope.c_str());
    auto results = store.GetGlobalResult( scope.c_str() );
    if(results.expired())
    {
        throw runtime_error("tplcombinationsolve: CombinePDF - Unable to retrieve tplcombination results");
    }
    std::shared_ptr<const CChisquareResult> result = std::dynamic_pointer_cast<const CChisquareResult>( results.lock() );

    CPdfz pdfz;
    Float64 cstLog = -1;

    Int32 retPdfz=-1;
    std::vector<TFloat64List> priors;
    std::vector<TFloat64List> chiSquares;
    std::vector<Float64> redshifts;
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
        if(cstLog==-1)
        {
            cstLog = result->CstLog;
            Log.LogInfo("tplcombinationsolve: using cstLog = %f", cstLog);
        }else if ( cstLog != result->CstLog)
        {
            Log.LogError("tplcombinationsolve: Found different cstLog values in results... val-1=%f != val-2=%f", cstLog, result->CstLog);
        }
        if(redshifts.size()==0)
        {
            redshifts = result->Redshifts;
        }

        for(Int32 kism=0; kism<nISM; kism++)
        {
            for(Int32 kigm=0; kigm<nIGM; kigm++)
            {
                priors.push_back(pdfz.GetConstantLogZPrior(result->Redshifts.size()));

                //correct chi2 for ampl. marg. if necessary: todo add switch, currently deactivated
                chiSquares.emplace_back(result->ChiSquareIntermediate.size(), DBL_MAX);
                for ( UInt32 kz=0; kz<result->Redshifts.size(); kz++)
                {
                    chiSquares.back()[kz] = result->ChiSquareIntermediate[kz][kism][kigm];// + resultXXX->ScaleMargCorrectionTplshapes[][]?;
                }
                Log.LogDetail("    tplcombinationsolve: Pdfz combine - prepared merit #%d for ism=%d, igm=%d", chiSquares.size()-1, kism, kigm);
            }
        }
    }

    if(chiSquares.size()>0)
    {
        if(opt_combine=="marg")
        {
            Log.LogInfo("    tplcombinationsolve: Pdfz combination - Marginalization");
            retPdfz = pdfz.Marginalize( redshifts, chiSquares, priors, cstLog, postmargZResult);
        }else if(opt_combine=="bestchi2")
        {
            Log.LogInfo("    tplcombinationsolve: Pdfz combination - BestChi2");
            retPdfz = pdfz.BestChi2( redshifts, chiSquares, priors, cstLog, postmargZResult);
        }else{
            Log.LogError("    tplcombinationsolve: Unable to parse pdf combination method option");
        }
    }else
    {
        Log.LogError("    tplcombinationsolve: Unable to find any chisquares prepared for combination. chiSquares.size()=%d", chiSquares.size());
    }


    if(retPdfz!=0)
    {
        Log.LogError("    tplcombinationsolve: Pdfz computation failed");
    }

    return retPdfz;
}

Bool CMethodTplcombinationSolve::ExtractCandidateResults(CDataStore &store, std::vector<Float64> zcandidates_unordered_list)
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



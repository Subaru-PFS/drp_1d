#include <RedshiftLibrary/method/tplcombinationsolve.h>
#include <RedshiftLibrary/method/templatefittingresult.h>
#include <RedshiftLibrary/operator/templatefittingresult.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/statistics/zprior.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <RedshiftLibrary/statistics/deltaz.h>
#include <RedshiftLibrary/spectrum/io/fitswriter.h>
#include <cfloat>

using namespace NSEpic;
using namespace std;


CMethodTplcombinationSolve::CMethodTplcombinationSolve()
{
}

CMethodTplcombinationSolve::~CMethodTplcombinationSolve()
{
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


std::shared_ptr<CTemplateFittingSolveResult> CMethodTplcombinationSolve::Compute(CDataStore& resultStore,
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
    m_radius = radius;
//    std::string _name = "Tplcombination";
    CDataStore::CAutoScope resultScope( resultStore, "tplcombinationsolve" );
//    std::string _scope = "tplcombination";
    std::string scopeStr = "chisquare";

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

    Solve(resultStore, spc, tplCatalog, tplCategoryList, lambdaRange, redshifts, overlapThreshold, maskList, _type, opt_interp, opt_extinction, opt_dustFit);

    storeResult = true;
//        }
//    }


    if( storeResult )
    {
        std::shared_ptr< CTemplateFittingSolveResult> solveResult = std::make_shared< CTemplateFittingSolveResult>(_type, "tplcombinationsolve");

        COperatorPdfz pdfz;

        {
            ChisquareArray chisquares = BuildChisquareArray(resultStore, scopeStr);
            pdfz.CombinePDF(chisquares, m_opt_pdfcombination);

            std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
            resultStore.StoreGlobalResult( pdfPath.c_str(), pdfz.m_postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
        }

        return solveResult;
    }

    return NULL;
}

Bool CMethodTplcombinationSolve::Solve(CDataStore& resultStore,
                                       const CSpectrum& spc,
                                       const CTemplateCatalog& tplCatalog,
                                       const TStringList& tplCategoryList,
                                       const TFloat64Range& lambdaRange,
                                       const TFloat64List& redshifts,
                                       Float64 overlapThreshold,
                                       std::vector<CMask> maskList,
                                       CTemplateFittingSolveResult::EType spctype,
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
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            tplList.push_back(tpl);
        }
    }


    //case: nType_all
    if(spctype == CTemplateFittingSolveResult::nType_all){
        _ntype = 3;
    }

    const CSpectrum::EType save_spcType = spc.GetType();

    std::vector<CSpectrum::EType> save_tplTypes;
    for ( CTemplate tpl: tplList){
        save_tplTypes.push_back(tpl.GetType());
    }

    for( Int32 i=0; i<_ntype; i++){
        if(spctype == CTemplateFittingSolveResult::nType_all){
            _spctype = _spctypetab[i];
        }else{
            _spctype = static_cast<CSpectrum::EType>(spctype);
        }

        spc.SetType(_spctype);
        for (CTemplate tpl: tplList)
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
            enable_dustFitting = 0;
        }

        // Compute merit function
        auto  result = std::dynamic_pointer_cast<CTemplateFittingResult>( m_tplcombinationOperator.Compute( spc, tplList, lambdaRange, redshifts, overlapThreshold, maskList, m_radius, opt_interp, enable_extinction, enable_dustFitting ) );

        if( !result )
        {
            //Log.LogError( "Failed to compute chi square value");
            return false;
        }else{
            // Store results
            Log.LogDetail("tplcombinationsolve: Save tplcombination results");
            resultStore.StoreScopedGlobalResult(scopeStr.c_str(), result );
            // Store spectrum results
            Log.LogDetail("tplcombinationsolve: Save spectrum/model results");
            m_tplcombinationOperator.SaveSpectrumResults(resultStore.getResultStore());
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

ChisquareArray CMethodTplcombinationSolve::BuildChisquareArray(const CDataStore& store, const std::string & scopeStr) const
{
    ChisquareArray chisquarearray;

    Log.LogDetail("tplcombinationsolve: build chisquare array");
    std::string scope = store.GetCurrentScopeName() + ".";
    scope.append(scopeStr.c_str());

    Log.LogDetail("    tplcombinationsolve: using results in scope: %s", scope.c_str());

    auto results = store.GetGlobalResult( scope.c_str() );
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

Bool CMethodTplcombinationSolve::ExtractCandidateResults(CDataStore& store, std::vector<Float64> zcandidates_unordered_list)
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

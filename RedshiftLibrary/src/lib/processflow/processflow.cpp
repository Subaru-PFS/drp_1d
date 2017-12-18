#include <RedshiftLibrary/processflow/processflow.h>

#include <RedshiftLibrary/continuum/standard.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/continuum/median.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/operator/peakdetection.h>
#include <RedshiftLibrary/gaussianfit/gaussianfit.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/operator/raymatching.h>
#include <RedshiftLibrary/common/median.h>

#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chicorr.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/method/blindsolve.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/operator/chisquare.h>

#include <RedshiftLibrary/operator/peakdetectionresult.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>

#include <RedshiftLibrary/method/chisquaresolve.h>
#include <RedshiftLibrary/method/chisquare2solve.h>
#include <RedshiftLibrary/method/chisquarelogsolve.h>
#include <RedshiftLibrary/method/correlationsolve.h>
#include <RedshiftLibrary/method/linematchingsolve.h>
#include <RedshiftLibrary/method/dtree7solve.h>
#include <RedshiftLibrary/method/dtree7solveresult.h>
#include <RedshiftLibrary/method/dtreebsolve.h>
#include <RedshiftLibrary/method/dtreebsolveresult.h>
#include <RedshiftLibrary/method/dtreecsolve.h>
#include <RedshiftLibrary/method/dtreecsolveresult.h>
#include <RedshiftLibrary/method/linematching2solve.h>
#include <RedshiftLibrary/method/linemodelsolve.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/method/linemodeltplshapesolve.h>
#include <RedshiftLibrary/method/linemodeltplshapesolveresult.h>

#include <RedshiftLibrary/reliability/zqual.h>

#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <float.h>

using namespace NSEpic;


CProcessFlow::CProcessFlow()
{

}

CProcessFlow::~CProcessFlow()
{

}

Bool CProcessFlow::Process( CProcessFlowContext& ctx )
{
    Log.LogInfo("<proc-spc><%s>", ctx.GetSpectrum().GetName().c_str());

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;

    ctx.GetParameterStore().Get( "lambdarange", lambdaRange );
    ctx.GetParameterStore().Get( "redshiftrange", redshiftRange );
    ctx.GetParameterStore().Get( "redshiftstep", redshiftStep );
    TFloat64Range spcLambdaRange;
    ctx.GetSpectrum().GetSpectralAxis().ClampLambdaRange( lambdaRange, spcLambdaRange );

    Log.LogInfo( "Processing spc:%s (LambdaRange: %f-%f:%f)", ctx.GetSpectrum().GetName().c_str(),
            spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
    DebugAssert( redshifts.size() > 0 );

    std::string CategoryFilter="all";
    // Remove Star category, and filter the list with regard to input variable CategoryFilter
    TStringList templateCategoryList;
    ctx.GetParameterStore().Get( "templateCategoryList", templateCategoryList );
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<templateCategoryList.size(); i++ )
    {
        std::string category = templateCategoryList[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    //retrieve the calibration dir path
    std::string calibrationDirPath;
    ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );


    Bool enableInputSpcCheck = true;
    if(enableInputSpcCheck)
    {
        //Check if the Spectrum is valid on the lambdarange
        const Float64 lmin = spcLambdaRange.GetBegin();
        const Float64 lmax = spcLambdaRange.GetEnd();
        if( !ctx.GetSpectrum().IsFluxValid( lmin, lmax ) ){
            Log.LogError( "Failed to validate spectrum flux: %s, on wavelength range (%.1f ; %.1f)", ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
            return false;
        }
        if( !ctx.GetSpectrum().IsNoiseValid( lmin, lmax ) ){
            Log.LogError( "Failed to validate noise from spectrum: %s, on wavelength range (%.1f ; %.1f)", ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
            return false;
        }
    }



    std::string methodName;
    ctx.GetParameterStore().Get( "method", methodName );
    boost::algorithm::to_lower(methodName);

    std::shared_ptr<COperatorResult> mResult;


    if(methodName  == "linemodel" ){

        CLineModelSolve Solve(calibrationDirPath);
        mResult = Solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 filteredTemplateCategoryList,
                                 ctx.GetRayCatalog(),
                                 spcLambdaRange,
                                 redshifts );



    }else if(methodName  == "chisquare2solve" ){
        Float64 overlapThreshold;
        ctx.GetParameterStore().Get( "chisquare2solve.overlapThreshold", overlapThreshold, 1.0);
        std::string opt_spcComponent;
        ctx.GetDataStore().GetScopedParam( "chisquare2solve.spectrum.component", opt_spcComponent, "raw" );
        std::string opt_interp;
        ctx.GetDataStore().GetScopedParam( "chisquare2solve.interpolation", opt_interp, "precomputedfinegrid" );
        std::string opt_extinction;
        ctx.GetDataStore().GetScopedParam( "chisquare2solve.extinction", opt_extinction, "no" );
        std::string opt_dustFit;
        ctx.GetDataStore().GetScopedParam( "chisquare2solve.dustfit", opt_dustFit, "no" );

        // prepare the unused masks
        std::vector<CMask> maskList;
        //retrieve the calibration dir path
        std::string calibrationDirPath;
        ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
        CMethodChisquare2Solve solve(calibrationDirPath);
        mResult = solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 filteredTemplateCategoryList,
                                 spcLambdaRange,
                                 redshifts,
                                 overlapThreshold,
                                 maskList,
                                 opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);

    }else if(methodName  == "chisquarelogsolve" ){
        Float64 overlapThreshold;
        ctx.GetParameterStore().Get( "chisquarelogsolve.overlapThreshold", overlapThreshold, 1.0);
        std::string opt_spcComponent;
        ctx.GetDataStore().GetScopedParam( "chisquarelogsolve.spectrum.component", opt_spcComponent, "raw" );
        std::string opt_interp="unused";
        std::string opt_extinction;
        ctx.GetDataStore().GetScopedParam( "chisquarelogsolve.extinction", opt_extinction, "no" );
        std::string opt_dustFit;
        ctx.GetDataStore().GetScopedParam( "chisquarelogsolve.dustfit", opt_dustFit, "no" );

        // prepare the unused masks
        std::vector<CMask> maskList;
        //retrieve the calibration dir path
        std::string calibrationDirPath;
        ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
        CMethodChisquareLogSolve solve(calibrationDirPath);
        mResult = solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 filteredTemplateCategoryList,
                                 spcLambdaRange,
                                 redshifts,
                                 overlapThreshold,
                                 maskList,
                                 opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);

    }else if(methodName  == "amazed0_1" ){
        COperatorDTree7Solve Solve(calibrationDirPath);
        mResult = Solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 templateCategoryList,
                                 ctx.GetRayCatalog(),
                                 lambdaRange,
                                 redshiftRange,
                                 redshiftStep);

    }else if(methodName  == "amazed0_2" ){
        COperatorDTreeBSolve Solve(calibrationDirPath);
        mResult = Solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 templateCategoryList,
                                 ctx.GetRayCatalog(),
                                 spcLambdaRange, redshifts);

    }else if(methodName  == "amazed0_3" ){
        COperatorDTreeCSolve Solve(calibrationDirPath);
        mResult = Solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 templateCategoryList,
                                 ctx.GetRayCatalog(),
                                 spcLambdaRange,
                                 redshifts);

    }else if(methodName  == "correlationsolve" ){
        COperatorCorrelationSolve solve;
        mResult = solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 filteredTemplateCategoryList,
                                 lambdaRange, redshiftRange, redshiftStep );

    }else if(methodName  == "blindsolve" ){
        COperatorBlindSolve blindSolve;
        mResult = blindSolve.Compute( ctx.GetDataStore(),
                                      ctx.GetSpectrum(),
                                      ctx.GetSpectrumWithoutContinuum(),
                                      ctx.GetTemplateCatalog(),
                                      filteredTemplateCategoryList,
                                      lambdaRange, redshiftRange, redshiftStep);

    }else if(methodName  == "linematching" ){
        COperatorLineMatchingSolve Solve;
        mResult = Solve.Compute(ctx.GetDataStore(), ctx.GetSpectrum(),
                                lambdaRange,
                                redshiftRange,
                                redshiftStep,
                                ctx.GetRayCatalog() );

    }else if(methodName  == "linematching2" ){
        COperatorLineMatching2Solve Solve;
        mResult = Solve.Compute(ctx.GetDataStore(),
                                ctx.GetSpectrum(),
                                spcLambdaRange,
                                redshiftRange,
                                redshiftStep,
                                ctx.GetRayCatalog() );

    }else{
        Log.LogError("Problem found while parsing the method parameter !");
        return false;
    }


    //Process Reliability estimation
    if(!mResult){
        Log.LogWarning( "Reliability skipped - no redshift results found");
    }else if(!isPdfValid(ctx)){
        Log.LogWarning( "Reliability skipped - no valid pdf result found");
    }else{
        CClassifierStore classifStore = ctx.GetClassifierStore();
        if(!classifStore.m_isInitialized)
        {
            Log.LogWarning( "Reliability not initialized. Skipped.");
        }else
        {
            Log.LogInfo( "Processing reliability");
            CQualz solve2;
            std::shared_ptr<const CQualzResult> solve2Result = solve2.Compute( ctx.GetDataStore(), classifStore, redshiftRange, redshiftStep );

            if(solve2Result)
            {
                std::string predLabel="";
                bool retPredLabel = solve2Result->GetPredictedLabel( ctx.GetDataStore(), predLabel );
                if( retPredLabel ) {
                    mResult->SetReliabilityLabel(predLabel);
                }else{
                    Log.LogError( "Unable estimate Reliability");
                }
            }else{
                Log.LogInfo( "No Reliability Result Found");
            }
        }
    }

    //finally save the method results with (optionally) the zqual label
    if( mResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", mResult );
    }else{
        Log.LogError( "Unable to store method result.");
        return false;
    }

    return true;
}

/**
 * @brief isPdfValid
 * @return
 */
Bool CProcessFlow::isPdfValid(CProcessFlowContext& ctx) const
{
    std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
    auto results_pdf =  ctx.GetDataStore().GetGlobalResult( scope_res.c_str() );
    auto logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results_pdf.lock() );

    if(!logzpdf1d)
    {
        return false;
    }

    if(logzpdf1d->Redshifts.size()<2)
    {
        return false;
    }

    //is it completely flat ?
    Float64 minVal=DBL_MAX;
    Float64 maxVal=-DBL_MAX;
    for(Int32 k=0; k<logzpdf1d->valProbaLog.size(); k++)
    {
        if(logzpdf1d->valProbaLog[k]<minVal)
        {
            minVal = logzpdf1d->valProbaLog[k];
        }
        if(logzpdf1d->valProbaLog[k]>maxVal)
        {
            maxVal = logzpdf1d->valProbaLog[k];
        }
    }
    if(minVal==maxVal){
        Log.LogError("PDF is flat !");
        return false;
    }

    return true;
}

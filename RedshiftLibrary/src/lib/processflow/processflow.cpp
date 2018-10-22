#include <RedshiftLibrary/processflow/processflow.h>

#include <RedshiftLibrary/continuum/standard.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>

#include <RedshiftLibrary/linemodel/calibrationconfig.h>

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
#include <RedshiftLibrary/method/zweimodelsolve.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/method/linemodeltplshapesolve.h>
#include <RedshiftLibrary/method/linemodeltplshapesolveresult.h>
#include <RedshiftLibrary/method/tplcombinationsolve.h>
#include <RedshiftLibrary/processflow/classificationresult.h>

#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/reliability/zqual.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <float.h>

using namespace std;
using namespace boost;
using namespace NSEpic;
namespace bfs = boost::filesystem;

CProcessFlow::CProcessFlow()
{

}

CProcessFlow::~CProcessFlow()
{

}

void CProcessFlow::Process( CProcessFlowContext& ctx )
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

    Log.LogInfo( "Processing spc:%s (CLambdaRange: %f-%f:%f)", ctx.GetSpectrum().GetName().c_str(),
            spcLambdaRange.GetBegin(), spcLambdaRange.GetEnd(), ctx.GetSpectrum().GetResolution());


    std::string methodName;
    ctx.GetParameterStore().Get( "method", methodName );
    boost::algorithm::to_lower(methodName);

    // Create redshift initial list by spanning redshift acdross the given range, with the given delta
    std::string redshiftSampling;
    ctx.GetParameterStore().Get( "redshiftsampling", redshiftSampling, "lin" ); //TODO: sampling in log cannot be used for now as zqual descriptors assume constant dz.
    TFloat64List raw_redshifts;
    if(redshiftSampling=="log")
    {
        raw_redshifts = redshiftRange.SpreadOverLog( redshiftStep ); //experimental: spreadover a grid at delta/(1+z), unusable because PDF needs regular z-step
    }else
    {
        raw_redshifts = redshiftRange.SpreadOver( redshiftStep );
    }


    //Override z-search grid for line measurement: load the zref values from a tsv catalog file (col0: spc name, col1: zref float value)
    std::string opt_linemeas_catalog_path;
    ctx.GetParameterStore().Get( "linemeascatalog", opt_linemeas_catalog_path, "" );
    TFloat64List redshifts;
    if(opt_linemeas_catalog_path!="")
    {
        Log.LogInfo( "Override z-search range !");
        Float64 zref = -1.0;
        namespace fs = boost::filesystem;
        Int32 reverseInclusionForIdMatching = 0; //0: because the names must match exactly, but: linemeas catalog includes the extension (.fits) and spc.GetName doesn't.
        bool computeOnZrange=false; //nb: hardcoded option for now
        Int32 colId = 2;//starts at 1, so that (for the linemeas_catalog) id_column=1, zref_column=2
        fs::path refFilePath(opt_linemeas_catalog_path.c_str());


        if ( fs::exists(refFilePath) )
        {
            std::string spcSubStringId = ctx.GetSpectrum().GetName();
            Log.LogInfo( "Override z-search: using spc-string: %s", spcSubStringId.c_str());
            getValueFromRefFile( refFilePath.c_str(), spcSubStringId, colId, zref, reverseInclusionForIdMatching);
        }
        if(zref==-1)
        {
          throw std::runtime_error("Override z-search: unable to find zref!");
        }

        if(computeOnZrange) //computing only on zref, or on a zrange around zref
        {
            Float64 deltaZrangeHalf = 0.5e-2; //override zrange
            Float64 stepZ = 1e-5;
            Float64 nStepsZ = deltaZrangeHalf*2/stepZ+1;
            for(Int32 kz=0; kz<nStepsZ; kz++)
            {
                Float64 _z = zref + kz*stepZ - deltaZrangeHalf;
                redshifts.push_back(_z);
            }
            Log.LogInfo( "Override z-search: zmin=%.5f, zmax=%.5f, zstep=%.5f", redshifts[0], redshifts[redshifts.size()-1], stepZ);
        }else{
            redshifts.push_back(zref);
        }

        if(methodName=="linemodel"){
            ctx.GetDataStore().SetScopedParam( "linemodelsolve.linemodel.extremacount", 1.0);
            Log.LogInfo( "Override z-search: Using overriden linemodelsolve.linemodel.extremacount: %f", 1.0);
            ctx.GetDataStore().SetScopedParam( "linemodelsolve.linemodel.firstpass.largegridstep", 0.0);
            Log.LogInfo( "Override z-search: Using overriden linemodelsolve.linemodel.firstpass.largegridstep: %f", 0.0);
        }

        Log.LogInfo( "Override z-search: Using overriden zref for spc %s : zref=%f", ctx.GetSpectrum().GetName().c_str(), zref);
    }else{
        redshifts = raw_redshifts;
    }

    if(redshifts.size() < 1)
    {
      Log.LogError("Unable to initialize the z-search grid (size=%d)", redshifts.size());
      throw std::runtime_error("Unable to initialize the z-search grid");
    }

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


    //************************************
    Bool enableInputSpcCorrect = true;
    std::string enableInputSpcCorrectStr;
    ctx.GetParameterStore().Get( "autocorrectinput", enableInputSpcCorrectStr, "no" );
    if(enableInputSpcCorrectStr!="yes")
    {
        enableInputSpcCorrect=false;
    }
    if(methodName  == "reliability" )
    {
        enableInputSpcCorrect = false;
    }
    if(enableInputSpcCorrect)
    {
        //Check if the Spectrum is valid on the lambdarange
        const Float64 lmin = spcLambdaRange.GetBegin();
        const Float64 lmax = spcLambdaRange.GetEnd();
        if( ctx.correctSpectrum( lmin, lmax ) ){
            Log.LogInfo( "Successfully corrected noise from spectrum: %s, on wavelength range (%.1f ; %.1f)", ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
        }
    }

    //************************************
    Bool enableInputSpcCheck = true;
    if(methodName  == "reliability" )
    {
        enableInputSpcCheck = false;
    }
    if(enableInputSpcCheck)
    {
        //Check if the Spectrum is valid on the lambdarange
        const Float64 lmin = spcLambdaRange.GetBegin();
        const Float64 lmax = spcLambdaRange.GetEnd();
        if( !ctx.GetSpectrum().IsFluxValid( lmin, lmax ) ){
            Log.LogError("Failed to validate spectrum flux: %s, on wavelength range (%.1f ; %.1f)",
                         ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
            throw std::runtime_error("Failed to validate spectrum flux");
        }else{
            Log.LogDetail( "Successfully validated spectrum flux: %s, on wavelength range (%.1f ; %.1f)", ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
        }
        if( !ctx.GetSpectrum().IsNoiseValid( lmin, lmax ) ){
            Log.LogError("Failed to validate noise from spectrum: %s, on wavelength range (%.1f ; %.1f)",
                         ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
            throw std::runtime_error("Failed to validate noise from spectrum");
        }else{
            Log.LogDetail( "Successfully validated noise from spectrum: %s, on wavelength range (%.1f ; %.1f)", ctx.GetSpectrum().GetName().c_str(), lmin, lmax );
        }
    }


    // Stellar method
    std::shared_ptr<COperatorResult> starResult;
    std::string enableStarFitting;
    ctx.GetParameterStore().Get( "enablestellarsolve", enableStarFitting, "no" );
    Log.LogInfo( "Stellar solve enabled : %s", enableStarFitting.c_str());
    if(enableStarFitting=="yes"){
        CDataStore::CAutoScope resultScope( ctx.GetDataStore(), "stellarsolve" );

        std::string calibrationDirPath;
        ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );

        bfs::path calibrationFolder( calibrationDirPath.c_str() );
        CCalibrationConfigHelper calibrationConfig;
        calibrationConfig.Init(calibrationDirPath);

	std::string starTemplates = calibrationConfig.Get_starTemplates_relpath();

            Log.LogInfo( "    Processflow - Loading star templates catalog : %s", starTemplates.c_str());
            std::string templateDir = (calibrationFolder/starTemplates.c_str()).string();

            TStringList   filteredStarTemplateCategoryList;
            filteredStarTemplateCategoryList.push_back( "star" );

            //temporary star catalog handling through calibration files, should be loaded somewhere else ?
            std::string medianRemovalMethod="zero";
            Float64 opt_medianKernelWidth = 150; //not used
            Int64 opt_nscales=8; //not used
            std::string dfBinPath="absolute_path_to_df_binaries_here"; //not used
            std::shared_ptr<CTemplateCatalog> starTemplateCatalog = std::shared_ptr<CTemplateCatalog>( new CTemplateCatalog(medianRemovalMethod, opt_medianKernelWidth, opt_nscales, dfBinPath) );
            starTemplateCatalog->Load( templateDir.c_str() );

	    for( UInt32 i=0; i<filteredStarTemplateCategoryList.size(); i++ )
	        {
                    std::string category = filteredStarTemplateCategoryList[i];
                    UInt32 ntpl = starTemplateCatalog->GetTemplateCount(category);
                    Log.LogInfo("Loaded (category=%s) template count = %d", category.c_str(), ntpl);
                }

            Float64 overlapThreshold;
            ctx.GetParameterStore().Get( "starsolve.overlapThreshold", overlapThreshold, 1.0);
            std::string opt_spcComponent;
            ctx.GetDataStore().GetScopedParam( "starsolve.spectrum.component", opt_spcComponent, "raw" );
            std::string opt_interp;
            ctx.GetDataStore().GetScopedParam( "starsolve.interpolation", opt_interp, "precomputedfinegrid" );
            std::string opt_extinction;
            ctx.GetDataStore().GetScopedParam( "starsolve.extinction", opt_extinction, "no" );
            std::string opt_dustFit;
            ctx.GetDataStore().GetScopedParam( "starsolve.dustfit", opt_dustFit, "yes" );

            // prepare the unused masks
            std::vector<CMask> maskList;
            //define the redshift search grid
            TFloat64Range starRedshiftRange=TFloat64Range(-0.5e-2, +0.5e-2);
            Float64 starRedshiftStep = 1e-5;
            Log.LogInfo("Stellar fitting redshift range = [%.5f, %.5f], step=%.6f", starRedshiftRange.GetBegin(), starRedshiftRange.GetEnd(), starRedshiftStep);
            TFloat64List stars_redshifts = starRedshiftRange.SpreadOver( starRedshiftStep );
            DebugAssert( stars_redshifts.size() > 0 );

            Log.LogInfo("Processing stellar fitting");
            //CMethodChisquare2Solve solve(calibrationDirPath);
            CMethodChisquareLogSolve solve(calibrationDirPath);
            starResult = solve.Compute( ctx.GetDataStore(),
                                     ctx.GetSpectrum(),
                                     ctx.GetSpectrumWithoutContinuum(),
                                     *starTemplateCatalog,
                                     filteredStarTemplateCategoryList,
                                     spcLambdaRange,
                                     stars_redshifts,
                                     overlapThreshold,
                                     maskList,
                                     "stellar_zPDF",
                                     opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);


            //finally save the stellar fitting results
            if( starResult ) {
                Log.LogInfo("Saving star fitting results");
                ctx.GetDataStore().StoreScopedGlobalResult( "stellarresult", starResult );
            }else{
                Log.LogError( "Unable to store stellar result.");
            }


    }

    // Galaxy method
    std::shared_ptr<COperatorResult> mResult;
    std::string galaxy_method_pdf_reldir = "zPDF";
    std::vector<Float64> zcandidates_unordered_list;
    if(methodName  == "linemodel" ){

        CLineModelSolve Solve(calibrationDirPath);
        mResult = Solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 filteredTemplateCategoryList,
                                 ctx.GetRayCatalog(),
                                 spcLambdaRange,
                                 redshifts,
                                 galaxy_method_pdf_reldir);

        if( mResult)
        {
            Log.LogInfo( "Extracting z-candidates from Linemodel method results" );
            std::shared_ptr<CLineModelSolveResult> solveResult = std::dynamic_pointer_cast<CLineModelSolveResult>( mResult );

            Bool retzc = solveResult->GetRedshiftCandidates( ctx.GetDataStore(), zcandidates_unordered_list);
            if(retzc)
            {
                Log.LogInfo( "Found %d z-candidates", zcandidates_unordered_list.size() );
            }else{
                Log.LogError( "Failed to get z candidates from these results");
            }
        }

    }else if(methodName  == "zweimodelsolve" ){

        CZweiModelSolve Solve(calibrationDirPath);
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
                                 galaxy_method_pdf_reldir,
                                 opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);

        if( mResult)
        {
            Log.LogInfo( "Extracting z-candidates from Chisquare2solve method results" );
            std::shared_ptr<CChisquare2SolveResult> solveResult = std::dynamic_pointer_cast<CChisquare2SolveResult>( mResult );
            Int32 n_cand = 5; //this is hardcoded for now for this method
            Bool retzc = solveResult->GetRedshiftCandidates( ctx.GetDataStore(), zcandidates_unordered_list, n_cand);
            if(retzc)
            {
                Log.LogInfo( "Found %d z-candidates", zcandidates_unordered_list.size() );
            }else{
                Log.LogError( "Failed to get z candidates from these results");
            }
        }

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
                                 galaxy_method_pdf_reldir,
                                 opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);


        if( mResult)
        {
            Log.LogInfo( "Extracting z-candidates from ChisquareLogsolve method results" );
            std::shared_ptr<CChisquareLogSolveResult> solveResult = std::dynamic_pointer_cast<CChisquareLogSolveResult>( mResult );
            Int32 n_cand = 5; //this is hardcoded for now for this method
            Bool retzc = solveResult->GetRedshiftCandidates( ctx.GetDataStore(), zcandidates_unordered_list, n_cand);
            if(retzc)
            {
                Log.LogInfo( "  Found %d z-candidates", zcandidates_unordered_list.size() );
            }else{
                Log.LogError( "  Failed to get z candidates from these results");
            }
        }


    }else if(methodName  == "tplcombinationsolve" ){
        Float64 overlapThreshold;
        ctx.GetParameterStore().Get( "tplcombinationsolve.overlapThreshold", overlapThreshold, 1.0);
        std::string opt_spcComponent;
        ctx.GetDataStore().GetScopedParam( "tplcombinationsolve.spectrum.component", opt_spcComponent, "raw" );
        std::string opt_interp="lin";
        ctx.GetDataStore().GetScopedParam( "tplcombinationsolve.interpolation", opt_interp, "lin" );
        std::string opt_extinction="no";
        //ctx.GetDataStore().GetScopedParam( "tplcombinationsolve.extinction", opt_extinction, "no" );
        std::string opt_dustFit="no";
        //ctx.GetDataStore().GetScopedParam( "tplcombinationsolve.dustfit", opt_dustFit, "no" );

        // prepare the unused masks
        std::vector<CMask> maskList;
        //retrieve the calibration dir path
        std::string calibrationDirPath;
        ctx.GetParameterStore().Get( "calibrationDir", calibrationDirPath );
        CMethodTplcombinationSolve solve(calibrationDirPath);
        mResult = solve.Compute( ctx.GetDataStore(),
                                 ctx.GetSpectrum(),
                                 ctx.GetSpectrumWithoutContinuum(),
                                 ctx.GetTemplateCatalog(),
                                 filteredTemplateCategoryList,
                                 spcLambdaRange,
                                 redshifts,
                                 overlapThreshold,
                                 maskList,
                                 galaxy_method_pdf_reldir,
                                 opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);

        if( mResult)
        {
            Log.LogInfo( "Extracting z-candidates from TplcombinationSolve method results" );
            std::shared_ptr<CTplcombinationSolveResult> solveResult = std::dynamic_pointer_cast<CTplcombinationSolveResult>( mResult );
            Int32 n_cand = 5; //this is hardcoded for now for this method
            Bool retzc = solveResult->GetRedshiftCandidates( ctx.GetDataStore(), zcandidates_unordered_list, n_cand);
            if(retzc)
            {
                Log.LogInfo( "Found %d z-candidates", zcandidates_unordered_list.size() );
            }else{
                Log.LogError( "Failed to get z candidates from these results");
            }
        }

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

    }else if(methodName  == "reliability" ){
        Log.LogInfo( "Processing RELIABILITY ONLY");
        //using an input pdf (ie. bypass redshift estimation method) from <intermSpcDir>/zPDF/logposterior.logMargP_Z_data.csv

        //loading pdf into datastore
        boost::filesystem::path perSpectrumDir="";
        perSpectrumDir = perSpectrumDir/( boost::filesystem::path( ctx.GetDataStore().GetProcessingID() ).string() );
        boost::filesystem::path inputPdfPath = perSpectrumDir/( boost::filesystem::path( "zPDF/logposterior.logMargP_Z_data.csv" ).string() ) ;

        Log.LogInfo( "Loading PDF from %s", inputPdfPath.string().c_str() );
        std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
        Int32 retPdfz = postmargZResult->Load(inputPdfPath.string().c_str());
        Log.LogInfo("Pdfz loaded n values = %d", postmargZResult->Redshifts.size());
        if(retPdfz<0)
        {
            Log.LogError("Pdfz loading failed (ret=%d)", retPdfz);
        }else{
            ctx.GetDataStore().StoreGlobalResult( "zPDF/logposterior.logMargP_Z_data", postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
        }
        mResult = std::shared_ptr<CLineModelSolveResult>(new CLineModelSolveResult());

    }else{
        throw std::runtime_error("Problem found while parsing the method parameter");
    }

    //Extracting the candidates summarized results (nb: only works for linemodel and chisquare methods as of Aug.2018)
    if(zcandidates_unordered_list.size()>0)
    {
        Log.LogInfo( "Computing candidates Probabilities" );
        std::shared_ptr<CPdfCandidateszResult> zcand = std::shared_ptr<CPdfCandidateszResult>(new CPdfCandidateszResult());

        std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
        auto results =  ctx.GetDataStore().GetGlobalResult( scope_res.c_str() );
        auto logzpdf1d = std::dynamic_pointer_cast<const CPdfMargZLogResult>( results.lock() );

        if(!logzpdf1d)
        {
            Log.LogError( "Extract Proba. for z candidates: no results retrieved from scope: %s", scope_res.c_str());
            throw std::runtime_error("Extract Proba. for z candidates: no results retrieved from scope");
        }

        Log.LogInfo( "  Integrating %d candidates proba.", zcandidates_unordered_list.size() );
        zcand->Compute(zcandidates_unordered_list, logzpdf1d->Redshifts, logzpdf1d->valProbaLog);
        ctx.GetDataStore().StoreScopedGlobalResult( "candidatesresult", zcand );
    }else{
        Log.LogInfo("No z-candidates found. Skipping probabilities integration...");
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
            CQualz solveReliab;
            std::shared_ptr<const CQualzResult> solveReliabResult = solveReliab.Compute( ctx.GetDataStore(), classifStore, redshiftRange, redshiftStep );

            if(solveReliabResult)
            {
                std::string predLabel="";
                bool retPredLabel = solveReliabResult->GetPredictedLabel( ctx.GetDataStore(), predLabel );
                if( retPredLabel ) {
                    mResult->SetReliabilityLabel(predLabel);
                }else{
                    Log.LogError( "Unable to estimate Reliability");
                }
            }else{
                Log.LogInfo( "No Reliability Result Found");
            }
        }
    }

    //estimate star/galaxy/qso classification
    Log.LogInfo("===============================================");
    std::shared_ptr<CClassificationResult> classifResult = std::shared_ptr<CClassificationResult>( new CClassificationResult() );
    Float64 qsoEvidence=-DBL_MAX;
    Float64 stellarEvidence=-DBL_MAX;
    Float64 galaxyEvidence=-DBL_MAX;
    mResult->GetEvidenceFromPdf(ctx.GetDataStore(), galaxyEvidence);
    Log.LogInfo( "Found galaxy evidence: %e", galaxyEvidence);
    std::string typeLabel = "G";
    if(enableStarFitting=="yes"){
        Int32 retStellarEv = starResult->GetEvidenceFromPdf(ctx.GetDataStore(), stellarEvidence);
        if(retStellarEv==0)
        {
            Log.LogInfo( "Found stellar evidence: %e", stellarEvidence);
            if(stellarEvidence>galaxyEvidence)
            {
                typeLabel = "S";
            }
        }
    }
    Log.LogInfo( "Setting object type: %s", typeLabel.c_str());
    mResult->SetTypeLabel(typeLabel); //maybe this is unecessary since there is a classifresult now
    classifResult->SetTypeLabel(typeLabel);
    classifResult->SetG(galaxyEvidence);
    classifResult->SetS(stellarEvidence);
    classifResult->SetQ(qsoEvidence);

    //save the classification results
    if( classifResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "classificationresult", classifResult );
    }else{
      throw std::runtime_error("Unable to store method result");
    }

    //finally save the method results with (optionally) the zqual label
    if( mResult ) {
        ctx.GetDataStore().StoreScopedGlobalResult( "redshiftresult", mResult );
    }else{
      throw std::runtime_error("Unable to store method result");
    }
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

    //is pdf any value nan ?
    for(Int32 k=0; k<logzpdf1d->valProbaLog.size(); k++)
    {
        if(logzpdf1d->valProbaLog[k] != logzpdf1d->valProbaLog[k])
        {
            Log.LogError("PDF has nan or invalid values !");
            return false;
        }
    }

    return true;
}


/**
 * \brief
 * Retrieve the true-redshift from a catalog file path
 *
 * reverseInclusion=0 (default): spcId is searched to be included in the Ref-File-Id
 * reverseInclusion=1 : Ref-File-Id is searched to be included in the spcId
 **/
Int32 CProcessFlow::getValueFromRefFile( const char* filePath, std::string spcid, Int32 colID, Float64& zref, Int32 reverseInclusion )
{
    ifstream file;

    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
        return false;

    string line;

    // Read file line by line
    while( getline( file, line ) )
    {
        // remove comments
        if(line.compare(0,1,"#",1)==0){
            continue;
        }
        char_separator<char> sep(" \t");

        // Tokenize each line
        typedef tokenizer< char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() && *it != "#" )
        {
            string name;
            if( it != tok.end() )
            {
                name = *it;
            }

            if(reverseInclusion==0)
            {
                std::size_t foundstr = name.find(spcid.c_str());
                if (foundstr==std::string::npos){
                    continue;
                }
            }else{
                std::size_t foundstr = spcid.find(name.c_str());
                if (foundstr==std::string::npos){
                    continue;
                }
            }

            // Found the correct spectrum ID: now read the ref values
            Int32 nskip = colID-1;
            for(Int32 i=0; i<nskip; i++)
            {
                ++it;
            }
            if( it != tok.end() )
            {

                zref = 0.0;
                try
                {
                    zref = lexical_cast<double>(*it);
                    return true;
                }
                catch (bad_lexical_cast&)
                {
                    zref = 0.0;
                    return false;
                }
            }

        }
    }
    file.close();
    return true;
}

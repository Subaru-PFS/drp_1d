#include <RedshiftLibrary/method/linemodelsolve.h>

#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/operator/pdfLogresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>


using namespace NSEpic;
using namespace std;
using namespace boost;


/** 
 * \brief Empty constructor.
 **/
CLineModelSolve::CLineModelSolve(string calibrationPath)
{
    m_calibrationPath = calibrationPath;
}

/**
 * \brief Empty destructor.
 **/
CLineModelSolve::~CLineModelSolve()
{

}

/**
 * \brief Returns a string describing the names and allowed values for the parameters of the Linemodel method.
 **/
const std::string CLineModelSolve::GetDescription()
{
    std::string desc;

    desc = "Method linemodel:\n";

    desc.append("\tparam: linemodel.linetypefilter = {""no"", ""E"", ""A""}\n");
    desc.append("\tparam: linemodel.lineforcefilter = {""no"", ""S""}\n");
    desc.append("\tparam: linemodel.fittingmethod = {""hybrid"", ""individual""}\n");
    desc.append("\tparam: linemodel.continuumcomponent = {""fromspectrum"", ""tplfit"", ""nocontinuum"", ""zero""}\n");
    desc.append("\tparam: linemodel.continuumismfit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.continuumigmfit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.continuumfitcount = <float value>\n");
    desc.append("\tparam: linemodel.secondpasslcfittingmethod = {""no"", ""svdlcp2""}\n");
    desc.append("\tparam: linemodel.rigidity = {""rules"", ""tplcorr"", ""tplshape""}\n");
    desc.append("\tparam: linemodel.tplratio_catalog = <relative path>\n");
    desc.append("\tparam: linemodel.tplratio_ismfit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.offsets_catalog = <relative path>\n");
    desc.append("\tparam: linemodel.linewidthtype = {""instrumentdriven"", ""velocitydriven"",  ""combined"",  ""nispvsspsf201707"", ""fixed""}\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.rules = {""all"", ""balmer"", ""strongweak"", ""superstrong"", ""ratiorange"", ""ciiiratio"", ""no""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");
    desc.append("\tparam: linemodel.extremacutprobathreshold = <float value> (-1:disabled)\n");
    desc.append("\tparam: linemodel.velocityfit = {""yes"", ""no""}\n");
    desc.append("\tparam: linemodel.emvelocityfitmin = <float value>\n");
    desc.append("\tparam: linemodel.emvelocityfitmax = <float value>\n");
    desc.append("\tparam: linemodel.emvelocityfitstep = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitmin = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitmax = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitstep = <float value>\n");
    //first pass
    desc.append("\tparam: linemodel.firstpass.largegridstep = <float value>, deactivated if negative or zero\n");
    desc.append("\tparam: linemodel.firstpass.tplratio_ismfit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.firstpass.multiplecontinuumfit_disable = {""no"", ""yes""}\n");

    desc.append("\tparam: linemodel.skipsecondpass = {""no"", ""yes""}\n");

    desc.append("\tparam: linemodel.pdfcombination = {""marg"", ""bestchi2""}\n");
    desc.append("\tparam: linemodel.stronglinesprior = <float value>, penalization factor = positive value or -1 to deactivate\n");
    desc.append("\tparam: linemodel.euclidnhaemittersStrength = <float value>, prior strength factor = positive value (typically 1 to 5) or -1 to deactivate\n");
    desc.append("\tparam: linemodel.saveintermediateresults = {""yes"", ""no""}\n");



    return desc;

}

/**
 * \brief
 * Populates the method parameters from the dataStore into the class members
 * Returns true if successful, false otherwise
 **/
Bool CLineModelSolve::PopulateParameters( CDataStore& dataStore )
{
    dataStore.GetScopedParam( "linemodel.linetypefilter", m_opt_linetypefilter, "no" );
    dataStore.GetScopedParam( "linemodel.lineforcefilter", m_opt_lineforcefilter, "no" );
    dataStore.GetScopedParam( "linemodel.fittingmethod", m_opt_fittingmethod, "hybrid" );
    dataStore.GetScopedParam( "linemodel.secondpasslcfittingmethod", m_opt_secondpasslcfittingmethod, "no" );
    dataStore.GetScopedParam( "linemodel.skipsecondpass", m_opt_skipsecondpass, "no" );
    dataStore.GetScopedParam( "linemodel.firstpass.fittingmethod", m_opt_firstpass_fittingmethod, "hybrid" );
    dataStore.GetScopedParam( "linemodel.firstpass.largegridstep", m_opt_firstpass_largegridstep, 0.001 );
    dataStore.GetScopedParam( "linemodel.firstpass.tplratio_ismfit", m_opt_firstpass_tplratio_ismfit, "no" );
    dataStore.GetScopedParam( "linemodel.firstpass.multiplecontinuumfit_disable", m_opt_firstpass_disablemultiplecontinuumfit, "yes" );

    std::string redshiftSampling;
    dataStore.GetParam( "redshiftsampling", redshiftSampling, "lin" ); //TODO: sampling in log cannot be used for now as zqual descriptors assume constant dz.
    if(redshiftSampling=="log")
    {
        m_opt_firstpass_largegridsampling = "log";
    }else{
        m_opt_firstpass_largegridsampling = "lin";
    }
    Log.LogDetail( "    firstpass - largegridsampling (auto set from redshiftsampling param.): %s", m_opt_firstpass_largegridsampling.c_str());

    dataStore.GetScopedParam( "linemodel.continuumcomponent", m_opt_continuumcomponent, "fromspectrum" );
    if(m_opt_continuumcomponent=="tplfit"){
        dataStore.GetScopedParam( "linemodel.continuumismfit", m_opt_tplfit_dustfit, "yes" );
        dataStore.GetScopedParam( "linemodel.continuumigmfit", m_opt_tplfit_igmfit, "yes" );
        dataStore.GetScopedParam( "linemodel.continuumfitcount", m_opt_continuumfitcount, 1 );
        dataStore.GetScopedParam( "linemodel.continuumfitignorelinesupport", m_opt_tplfit_ignoreLinesSupport, "no" );
    }
    dataStore.GetScopedParam( "linemodel.rigidity", m_opt_rigidity, "rules" );
    if(m_opt_rigidity=="tplshape")
    {
        dataStore.GetScopedParam( "linemodel.tplratio_catalog", m_opt_tplratio_reldirpath, "linecatalogs_tplshapes/linecatalogs_tplshape_ExtendedTemplatesJan2017v3_20170602_B14C_v3_emission" );
        dataStore.GetScopedParam( "linemodel.tplratio_ismfit", m_opt_tplratio_ismfit, "yes" );
    }
    dataStore.GetScopedParam( "linemodel.offsets_catalog", m_opt_offsets_reldirpath, "linecatalogs_offsets/offsetsCatalogs_20170410_m150" );

    dataStore.GetScopedParam( "linemodel.linewidthtype", m_opt_lineWidthType, "velocitydriven" );
    dataStore.GetScopedParam( "linemodel.instrumentresolution", m_opt_resolution, 2350.0 );
    dataStore.GetScopedParam( "linemodel.velocityemission", m_opt_velocity_emission, 100.0 );
    dataStore.GetScopedParam( "linemodel.velocityabsorption", m_opt_velocity_absorption, 300.0 );
    dataStore.GetScopedParam( "linemodel.velocityfit", m_opt_velocityfit, "yes" );
    if(m_opt_velocityfit=="yes"){
        dataStore.GetScopedParam( "linemodel.emvelocityfitmin", m_opt_em_velocity_fit_min, 20.0 );
        dataStore.GetScopedParam( "linemodel.emvelocityfitmax", m_opt_em_velocity_fit_max, 500.0 );
        dataStore.GetScopedParam( "linemodel.emvelocityfitstep", m_opt_em_velocity_fit_step, 20.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitmin", m_opt_abs_velocity_fit_min, 150.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitmax", m_opt_abs_velocity_fit_max, 500.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitstep", m_opt_abs_velocity_fit_step, 20.0 );
    }
    dataStore.GetScopedParam( "linemodel.continuumreestimation", m_opt_continuumreest, "no" );
    dataStore.GetScopedParam( "linemodel.rules", m_opt_rules, "all" );
    dataStore.GetScopedParam( "linemodel.extremacount", m_opt_extremacount, 10.0 );
    dataStore.GetScopedParam( "linemodel.extremacutprobathreshold", m_opt_candidatesLogprobaCutThreshold, -1 );
    dataStore.GetScopedParam( "linemodel.stronglinesprior", m_opt_stronglinesprior, -1);
    dataStore.GetScopedParam( "linemodel.euclidnhaemittersStrength", m_opt_euclidNHaEmittersPriorStrength, -1);
    dataStore.GetScopedParam( "linemodel.pdfcombination", m_opt_pdfcombination, "marg");
    dataStore.GetScopedParam( "linemodel.saveintermediateresults", m_opt_saveintermediateresults, "no");

    //Auto-correct fitting method
    std::string forcefittingmethod = "individual";
    if(m_opt_rigidity=="tplshape" && m_opt_fittingmethod == "hybrid")
    {
        m_opt_fittingmethod = forcefittingmethod;
        dataStore.SetScopedParam("linemodel.fittingmethod", m_opt_fittingmethod);
        Log.LogInfo( "LineModel fitting method auto-correct due to tplshape rigidity");
    }
    if(m_opt_rigidity=="tplshape" && m_opt_firstpass_fittingmethod == "hybrid")
    {
        m_opt_firstpass_fittingmethod = forcefittingmethod;
        dataStore.SetScopedParam("linemodel.firstpass.fittingmethod", m_opt_firstpass_fittingmethod);
        Log.LogInfo( "LineModel first pass fitting method auto-correct due to tplshape rigidity");
    }
    Log.LogInfo( "Linemodel parameters:");
    Log.LogInfo( "    -linetypefilter: %s", m_opt_linetypefilter.c_str());
    Log.LogInfo( "    -lineforcefilter: %s", m_opt_lineforcefilter.c_str());
    Log.LogInfo( "    -fittingmethod: %s", m_opt_fittingmethod.c_str());
    Log.LogInfo( "    -linewidthtype: %s", m_opt_lineWidthType.c_str());
    if(m_opt_lineWidthType=="combined"){
        Log.LogInfo( "    -instrumentresolution: %.2f", m_opt_resolution);
        Log.LogInfo( "    -velocity emission: %.2f", m_opt_velocity_emission);
        Log.LogInfo( "    -velocity absorption: %.2f", m_opt_velocity_absorption);
        Log.LogInfo( "    -velocity fit: %s", m_opt_velocityfit.c_str());
    }else if(m_opt_lineWidthType=="instrumentdriven"){
        Log.LogInfo( "    -instrumentresolution: %.2f", m_opt_resolution);
    }else if(m_opt_lineWidthType=="velocitydriven"){
        Log.LogInfo( "    -velocity emission: %.2f", m_opt_velocity_emission);
        Log.LogInfo( "    -velocity absorption: %.2f", m_opt_velocity_absorption);
        Log.LogInfo( "    -velocity fit: %s", m_opt_velocityfit.c_str());
    }else if(m_opt_lineWidthType=="nispvsspsf201707"){
        Log.LogInfo( "    -source size: hardcoded");
        Log.LogInfo( "    -velocity emission: %.2f", m_opt_velocity_emission);
        Log.LogInfo( "    -velocity absorption: %.2f", m_opt_velocity_absorption);
        Log.LogInfo( "    -velocity fit: %s", m_opt_velocityfit.c_str());
    }else if(m_opt_lineWidthType=="nispsim2016"){
        Log.LogInfo( "    -velocity emission: %.2f", m_opt_velocity_emission);
        Log.LogInfo( "    -velocity absorption: %.2f", m_opt_velocity_absorption);
        Log.LogInfo( "    -velocity fit: %s", m_opt_velocityfit.c_str());
    }
    if(m_opt_velocityfit=="yes"){
        Log.LogInfo( "    -em velocity fit min : %.1f", m_opt_em_velocity_fit_min);
        Log.LogInfo( "    -em velocity fit max : %.1f", m_opt_em_velocity_fit_max);
        Log.LogInfo( "    -em velocity fit step : %.1f", m_opt_em_velocity_fit_step);
        Log.LogInfo( "    -abs velocity fit min : %.1f", m_opt_abs_velocity_fit_min);
        Log.LogInfo( "    -abs velocity fit max : %.1f", m_opt_abs_velocity_fit_max);
        Log.LogInfo( "    -abs velocity fit step : %.1f", m_opt_abs_velocity_fit_step);
    }

    Log.LogInfo( "    -rigidity: %s", m_opt_rigidity.c_str());
    if(m_opt_rigidity=="rules"){
        Log.LogInfo( "      -rules: %s", m_opt_rules.c_str());
    }else if(m_opt_rigidity=="tplshape")
    {
        Log.LogInfo( "      -tplratio_catalog: %s", m_opt_tplratio_reldirpath.c_str());
        Log.LogInfo( "      -tplratio_ismfit: %s", m_opt_tplratio_ismfit.c_str());
    }
    Log.LogInfo( "    -linemodel offsets_catalog: %s", m_opt_offsets_reldirpath.c_str());

    Log.LogInfo( "    -continuumcomponent: %s", m_opt_continuumcomponent.c_str());
    if(m_opt_continuumcomponent=="tplfit"){
        Log.LogInfo( "      -tplfit_ismfit: %s", m_opt_tplfit_dustfit.c_str());
        Log.LogInfo( "      -tplfit_igmfit: %s", m_opt_tplfit_igmfit.c_str());
        Log.LogInfo( "      -continuum fit count:  %.0f", m_opt_continuumfitcount);
        Log.LogInfo( "      -tplfit_ignorelinesupport: %s", m_opt_tplfit_ignoreLinesSupport.c_str());
        Log.LogInfo( "      -tplfit_secondpass-LC-fitting-method: %s", m_opt_secondpasslcfittingmethod.c_str());
    }
    Log.LogInfo( "    -continuumreestimation: %s", m_opt_continuumreest.c_str());
    Log.LogInfo( "    -extremacount: %.0f", m_opt_extremacount);
    Log.LogInfo( "    -extrema cut proba-threshold: %.0f", m_opt_candidatesLogprobaCutThreshold);
    Log.LogInfo( "    -first pass:");
    Log.LogInfo( "      -largegridstep: %.6f", m_opt_firstpass_largegridstep);
    Log.LogInfo( "      -fittingmethod: %s", m_opt_firstpass_fittingmethod.c_str());
    Log.LogInfo( "      -tplratio_ismfit: %s", m_opt_firstpass_tplratio_ismfit.c_str());
    Log.LogInfo( "      -multiplecontinuumfit_disable: %s", m_opt_firstpass_disablemultiplecontinuumfit.c_str());


    Log.LogInfo( "    -skip second pass: %s", m_opt_skipsecondpass.c_str());

    Log.LogInfo( "    -pdf-stronglinesprior: %e", m_opt_stronglinesprior);
    Log.LogInfo( "    -pdf-euclidNHaEmittersPriorStrength: %e", m_opt_euclidNHaEmittersPriorStrength);
    Log.LogInfo( "    -pdf-combination: %s", m_opt_pdfcombination.c_str()); // "marg";    // "bestchi2";    // "bestproba";

    if(m_opt_saveintermediateresults=="yes")
    {
        m_opt_enableSaveChisquareTplshapeResults = true;
    }else{
        m_opt_enableSaveChisquareTplshapeResults = false;
    }
    Log.LogInfo( "    -save-intermediate-chisquaretplshaperesults: %d", (int)m_opt_enableSaveChisquareTplshapeResults);


    return true;
}

/**
 * \brief Calls the Solve method and returns a new "result" object.
 * Call Solve.
 * Return a pointer to an empty CLineModelSolveResult. (The results for Linemodel will reside in the linemodel.linemodel result).
 **/
std::shared_ptr<CLineModelSolveResult> CLineModelSolve::Compute( CDataStore& dataStore,
                                                                       const CSpectrum& spc,
                                                                       const CSpectrum& spcWithoutCont,
                                       const CTemplateCatalog& tplCatalog,
                                       const TStringList& tplCategoryList,
                                                                       const CRayCatalog& restraycatalog,
                                                                       const TFloat64Range& lambdaRange,
                                       const TFloat64List& redshifts,
                                       const std::string outputPdfRelDir)
{
    CDataStore::CAutoScope resultScope( dataStore, "linemodelsolve" );

    PopulateParameters( dataStore );
    Int32 retSolve = Solve( dataStore, spc, spcWithoutCont, tplCatalog, tplCategoryList, restraycatalog, lambdaRange, redshifts);

    if(retSolve){

        /* ------------------------  COMPUTE POSTERIOR PDF  --------------------------  */
        Log.LogInfo("linemodelsolve: Pdfz computation");

        std::string scope = "linemodelsolve.linemodel";
        auto results = dataStore.GetGlobalResult( scope.c_str() );
        if(results.expired())
        {
            Log.LogError("linemodelsolve: Unable to retrieve linemodel results");
            return NULL;
        }
        std::shared_ptr<const CLineModelResult> result = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );

        std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
        std::shared_ptr<CPdfLogResult> zpriorResult = std::shared_ptr<CPdfLogResult>(new CPdfLogResult());
        Int32 retCombinePdf = CombinePDF(result, m_opt_rigidity, m_opt_pdfcombination, m_opt_stronglinesprior, m_opt_euclidNHaEmittersPriorStrength, postmargZResult, zpriorResult);

        if(retCombinePdf!=0)
        {
            Log.LogError("Linemodel: Pdfz computation failed");
        }else{

            Log.LogDetail("    linemodelsolve: Storing priors (size=%d)", zpriorResult->Redshifts.size());
            std::string priorPath = outputPdfRelDir+"/logprior.logP_Z_data";
            dataStore.StoreGlobalResult( priorPath.c_str(), zpriorResult);

            //check pdf sum=1
            CPdfz pdfz;
            Float64 sumRect = pdfz.getSumRect(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Float64 sumTrapez = pdfz.getSumTrapez(postmargZResult->Redshifts, postmargZResult->valProbaLog);
            Log.LogDetail("    linemodelsolve: Pdfz normalization - sum rect. = %e", sumRect);
            Log.LogDetail("    linemodelsolve: Pdfz normalization - sum trapz. = %e", sumTrapez);
            Bool pdfSumCheck = abs(sumRect-1.0)<1e-1 || abs(sumTrapez-1.0)<1e-1;
            if(!pdfSumCheck){
                Log.LogWarning("    linemodelsolve: Pdfz normalization failed (rectsum = %f, trapzesum = %f)", sumRect, sumTrapez);
            }

            Log.LogDetail("    linemodelsolve: Storing PDF results");
            std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
            dataStore.StoreGlobalResult( pdfPath.c_str(), postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
        }


        //SaveContinuumPDF(dataStore, result);

        //Save chisquareTplshape results
        if(m_opt_enableSaveChisquareTplshapeResults)
        {
            for(Int32 km=0; km<result->ChiSquareTplshapes.size(); km++)
            {
                std::shared_ptr<CLineModelResult> result_chisquaretplshape = std::shared_ptr<CLineModelResult>( new CLineModelResult() );
                result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0, std::vector<Float64>());
                for(Int32 kz=0; kz<result->Redshifts.size(); kz++)
                {
                    result_chisquaretplshape->ChiSquare[kz] = result->ChiSquareTplshapes[km][kz];
                }

                std::string resname = (boost::format("linemodel_chisquaretplshape/linemodel_chisquaretplshape_%d") % km).str();
                dataStore.StoreScopedGlobalResult( resname.c_str(), result_chisquaretplshape );
            }


            //Save scaleMargCorrTplshape results
            for(Int32 km=0; km<result->ScaleMargCorrectionTplshapes.size(); km++)
            {
                std::shared_ptr<CLineModelResult> result_chisquaretplshape = std::shared_ptr<CLineModelResult>( new CLineModelResult() );
                result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0, std::vector<Float64>());
                for(Int32 kz=0; kz<result->Redshifts.size(); kz++)
                {
                    result_chisquaretplshape->ChiSquare[kz] = result->ScaleMargCorrectionTplshapes[km][kz];
                }

                std::string resname = (boost::format("linemodel_chisquaretplshape/linemodel_scalemargcorrtplshape_%d") % km).str();
                dataStore.StoreScopedGlobalResult( resname.c_str(), result_chisquaretplshape );
            }
        }
    }else{
        return NULL;
    }

    return std::shared_ptr<CLineModelSolveResult>( new CLineModelSolveResult() );
}


Int32 CLineModelSolve::CombinePDF(std::shared_ptr<const CLineModelResult> result,
                                  std::string opt_rigidity,
                                  std::string opt_combine,
                                  Float64 opt_stronglinesprior,
                                  Float64 opt_euclidNHaEmittersPriorStrength,
                                  std::shared_ptr<CPdfMargZLogResult> postmargZResult,
                                  std::shared_ptr<CPdfLogResult> zPrior)
{
    bool zPriorStrongLinePresence = (opt_stronglinesprior>0.0);
    if(zPriorStrongLinePresence)
    {
        Log.LogInfo("Linemodel: Pdfz computation: StrongLinePresence prior enabled: factor=%e", opt_stronglinesprior);
    }else{
        Log.LogInfo("Linemodel: Pdfz computation: StrongLinePresence prior disabled");
    }
    Float64 opt_nlines_snr_penalization_factor = -1;
    bool zPriorNLineSNR = (opt_nlines_snr_penalization_factor>0.0);
    if(zPriorNLineSNR)
    {
        Log.LogInfo("Linemodel: Pdfz computation: N lines snr>cut prior enabled: factor=%e", opt_nlines_snr_penalization_factor);
    }else{
        Log.LogInfo("Linemodel: Pdfz computation: N lines snr>cut prior disabled");
    }

    //hardcoded Euclid-NHaZprior parameter
    bool zPriorEuclidNHa = false;
    if(opt_euclidNHaEmittersPriorStrength>0.0)
    {
        zPriorEuclidNHa = true;
        Log.LogInfo("Linemodel: Pdfz computation: EuclidNHa prior enabled, with strength-coeff: %e", opt_euclidNHaEmittersPriorStrength);
    }else{
        Log.LogInfo("Linemodel: Pdfz computation: EuclidNHa prior disabled");
    }

    Log.LogInfo("Linemodel: Pdfz computation");
    CPdfz pdfz;
    Float64 cstLog = result->cstLog;
    TFloat64List logProba;
    Float64 logEvidence;

    Int32 retPdfz=-1;
    if(opt_rigidity!="tplshape" || opt_combine=="bestchi2")
    {
        if(opt_rigidity!="tplshape")
        {
            Log.LogInfo("Linemodel: Pdfz computation - simple (no combination)");
        }
        if(opt_combine=="bestchi2")
        {
            Log.LogInfo("Linemodel: Pdfz computation - simple (method=bestchi2)");
        }
        zPrior->SetSize(result->Redshifts.size());
        for ( UInt32 k=0; k<result->Redshifts.size(); k++)
        {
            zPrior->Redshifts[k] = result->Redshifts[k];
        }
        if(zPriorStrongLinePresence)
        {
            UInt32 lineTypeFilter = 1;// for emission lines only
            std::vector<bool> strongLinePresence = result->GetStrongLinesPresence(lineTypeFilter, result->LineModelSolutions);
            //std::vector<bool> strongLinePresence = result->GetStrongestLineIsHa(result->LineModelSolutions); //warning: hardcoded selpp replaced by whasp for lm-tplratio
            zPrior->valProbaLog = pdfz.GetStrongLinePresenceLogZPrior(strongLinePresence, opt_stronglinesprior);
        }else{
            zPrior->valProbaLog = pdfz.GetConstantLogZPrior(result->Redshifts.size());
        }
        if(zPriorEuclidNHa)
        {
            std::vector<Float64> zlogPriorNHa = pdfz.GetEuclidNhaLogZPrior(result->Redshifts, opt_euclidNHaEmittersPriorStrength);
            zPrior->valProbaLog = pdfz.CombineLogZPrior(zPrior->valProbaLog, zlogPriorNHa);
        }
        if(zPriorNLineSNR)
        {
            std::vector<Int32> n_lines_above_snr = result->GetNLinesAboveSnrcut(result->LineModelSolutions);
            std::vector<Float64> zlogPriorNLinesAboveSNR = pdfz.GetNLinesSNRAboveCutLogZPrior(n_lines_above_snr, opt_nlines_snr_penalization_factor);
            zPrior->valProbaLog = pdfz.CombineLogZPrior(zPrior->valProbaLog, zlogPriorNLinesAboveSNR);
        }

        //correct chi2 if necessary: todo add switch
        TFloat64List logLikelihoodCorrected(result->ChiSquare.size(), DBL_MAX);
        for ( UInt32 k=0; k<result->Redshifts.size(); k++)
        {
            logLikelihoodCorrected[k] = result->ChiSquare[k];// + result->ScaleMargCorrection[k];
        }
        retPdfz = pdfz.Compute(logLikelihoodCorrected, result->Redshifts, cstLog, zPrior->valProbaLog, logProba, logEvidence);
        if(retPdfz==0){
            postmargZResult->countTPL = result->Redshifts.size(); // assumed 1 model per z
            postmargZResult->Redshifts.resize(result->Redshifts.size());
            postmargZResult->valProbaLog.resize(result->Redshifts.size());
            for ( UInt32 k=0; k<result->Redshifts.size(); k++)
            {
                postmargZResult->Redshifts[k] = result->Redshifts[k] ;
                postmargZResult->valProbaLog[k] = logProba[k];
            }
            postmargZResult->valEvidenceLog = logEvidence;
        }
    }
    else if(opt_combine=="bestproba" || opt_combine=="marg"){

        Log.LogInfo("Linemodel: Pdfz computation - combination: method=%s, n=%d", opt_combine.c_str(), result->ChiSquareTplshapes.size());
        std::vector<TFloat64List> zpriorsTplshapes;
        for(Int32 k=0; k<result->ChiSquareTplshapes.size(); k++)
        {
            TFloat64List _prior;
            if(zPriorStrongLinePresence)
            {
                std::vector<bool> strongLinePresence = result->StrongELPresentTplshapes[k];
                _prior = pdfz.GetStrongLinePresenceLogZPrior(strongLinePresence, opt_stronglinesprior);
            }else
            {
                _prior = pdfz.GetConstantLogZPrior(result->Redshifts.size());
            }

            if(zPriorEuclidNHa)
            {
                std::vector<Float64> zlogPriorNHa = pdfz.GetEuclidNhaLogZPrior(result->Redshifts, opt_euclidNHaEmittersPriorStrength);
                _prior = pdfz.CombineLogZPrior(_prior, zlogPriorNHa);
            }
            if(zPriorNLineSNR)
            {
                std::vector<Int32> n_lines_above_snr = result->NLinesAboveSNRTplshapes[k];
                std::vector<Float64> zlogPriorNLinesAboveSNR = pdfz.GetNLinesSNRAboveCutLogZPrior(n_lines_above_snr, opt_nlines_snr_penalization_factor);
                _prior = pdfz.CombineLogZPrior(_prior, zlogPriorNLinesAboveSNR);
            }
            zpriorsTplshapes.push_back(_prior);
        }

        //correct chi2 if necessary: todo add switch
        std::vector<TFloat64List> ChiSquareTplshapesCorrected;
        for(Int32 k=0; k<result->ChiSquareTplshapes.size(); k++)
        {
            TFloat64List logLikelihoodCorrected(result->ChiSquareTplshapes[k].size(), DBL_MAX);
            for ( UInt32 kz=0; kz<result->Redshifts.size(); kz++)
            {
                logLikelihoodCorrected[kz] = result->ChiSquareTplshapes[k][kz];// + result->ScaleMargCorrectionTplshapes[k][kz];
            }
            ChiSquareTplshapesCorrected.push_back(logLikelihoodCorrected);
        }

        if(opt_combine=="marg")
        {
            retPdfz = pdfz.Marginalize( result->Redshifts, ChiSquareTplshapesCorrected, zpriorsTplshapes, cstLog, postmargZResult, result->PriorTplshapes);
        }else{
            retPdfz = pdfz.BestProba( result->Redshifts, ChiSquareTplshapesCorrected, zpriorsTplshapes, cstLog, postmargZResult);
        }
        // todo: store priors for each tplshape model ?
    }else{
        Log.LogError("Linemodel: Unable to parse pdf combination method option");
    }

    return retPdfz;
}


Int32 CLineModelSolve::SaveContinuumPDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result)
{
    Log.LogInfo("Linemodel: continuum Pdfz computation");
    std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
    CPdfz pdfz;
    Float64 cstLog = result->cstLog;
    TFloat64List logProba;
    Float64 logEvidence;

    Int32 retPdfz=-1;

    std::shared_ptr<CPdfLogResult> zPrior = std::shared_ptr<CPdfLogResult>(new CPdfLogResult());
    zPrior->SetSize(result->Redshifts.size());
    for ( UInt32 k=0; k<result->Redshifts.size(); k++)
    {
        zPrior->Redshifts[k] = result->Redshifts[k];
    }

    zPrior->valProbaLog = pdfz.GetConstantLogZPrior(result->Redshifts.size());


    //correct chi2 if necessary: todo add switch
    TFloat64List logLikelihoodCorrected(result->ChiSquareContinuum.size(), DBL_MAX);
    for ( UInt32 k=0; k<result->Redshifts.size(); k++)
    {
        logLikelihoodCorrected[k] = result->ChiSquareContinuum[k];// + result->ScaleMargCorrectionContinuum[k];
    }
    retPdfz = pdfz.Compute(logLikelihoodCorrected, result->Redshifts, cstLog, zPrior->valProbaLog, logProba, logEvidence);
    if(retPdfz==0){
        store.StoreGlobalResult( "zPDF/logpriorcontinuum.logP_Z_data", zPrior);

        postmargZResult->countTPL = result->Redshifts.size(); // assumed 1 model per z
        postmargZResult->Redshifts.resize(result->Redshifts.size());
        postmargZResult->valProbaLog.resize(result->Redshifts.size());
        for ( UInt32 k=0; k<result->Redshifts.size(); k++)
        {
            postmargZResult->Redshifts[k] = result->Redshifts[k] ;
            postmargZResult->valProbaLog[k] = logProba[k];
        }

        store.StoreGlobalResult( "zPDF/logposteriorcontinuum.logMargP_Z_data", postmargZResult);
    }else{
        Log.LogError("Linemodel: Pdfz computation for continuum FAILED!");
    }

    return 0;
}



/**
 * \brief
 * Retrieve the true-velocities from a hardcoded ref file path
 * nb: this is a hack for development purposes
 **/
Int32 getVelocitiesFromRefFile( const char* filePath, std::string spcid, Float64& elv, Float64& alv )
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
            std::size_t foundstr = name.find(spcid.c_str());
            if (foundstr==std::string::npos){
                continue;
            }

            // Found the correct spectrum ID: now read the ref values
            Int32 nskip = 7;
            for(Int32 i=0; i<nskip; i++)
            {
                ++it;
            }
            if( it != tok.end() )
            {

                elv = 0.0;
                try
                {
                    elv = lexical_cast<double>(*it);
                }
                catch (bad_lexical_cast&)
                {
                    elv = 0.0;
                    return false;
                }
            }
            ++it;
            if( it != tok.end() )
            {
                alv = 0.0;
                try
                {
                    alv = lexical_cast<double>(*it);
                }
                catch (bad_lexical_cast&)
                {
                    alv = 0.0;
                    return false;
                }

            }
        }
    }
    file.close();
    return true;
}

/**
 * \brief
 * Create a continuum object by subtracting spcWithoutCont from the spc.
 * Configure the opt_XXX variables from the dataStore scope parameters.
 * LogInfo the opt_XXX values.
 * Create a COperatorLineModel, call its Compute method. 
 * If that returned true, store results.
 **/
Bool CLineModelSolve::Solve( CDataStore& dataStore,
                             const CSpectrum& spc,
                             const CSpectrum& spcWithoutCont,
                             const CTemplateCatalog& tplCatalog,
                             const TStringList& tplCategoryList,
                             const CRayCatalog& restraycatalog,
                             const TFloat64Range& lambdaRange,
                             const TFloat64List& redshifts )
{
    std::string scopeStr = "linemodel";

    CSpectrum _spc = spc;
    CSpectrum _spcContinuum = spc;
    _spcContinuum.SetMedianWinsize(spcWithoutCont.GetMedianWinsize());
    _spcContinuum.SetDecompScales(spcWithoutCont.GetDecompScales());
    _spcContinuum.SetContinuumEstimationMethod(spcWithoutCont.GetContinuumEstimationMethod());
    _spcContinuum.SetWaveletsDFBinPath(spcWithoutCont.GetWaveletsDFBinPath());
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract( spcWithoutCont.GetFluxAxis() );
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;

    //    //Hack: load the simulated true-velocities
    //    if(false)
    //    {
    //        Float64 elv = 0.0;
    //        Float64 alv = 0.0;
    //        namespace fs = boost::filesystem;
    //        fs::path refFilePath("/home/aschmitt/data/simu_linemodel/simulm_20160513/simulation_pfswlinemodel_20160513_10spcperbin/refz.txt");
    //        if ( fs::exists(refFilePath) )
    //        {
    //            std::string spcSubStringId = spc.GetName().substr(0, 20);
    //            getVelocitiesFromRefFile( refFilePath.c_str(), spcSubStringId, elv, alv);
    //        }
    //        Float64 offsetv = 0.0;
    //        m_opt_velocity_emission = elv+offsetv;
    //        m_opt_velocity_absorption = alv+offsetv;
    //        Log.LogInfo( "Linemodel - hack - Loaded velocities for spc %s : elv=%4.1f, alv=%4.1f", spc.GetName().c_str(), elv, alv);
    //    }

    // Compute with linemodel operator
    COperatorLineModel linemodel;
    Int32 retInit = linemodel.Init(_spc, redshifts);
    if( retInit!=0 )
    {
        Log.LogError( "Line Model, init failed. Aborting" );
        return false;
    }
    linemodel.m_opt_firstpass_fittingmethod=m_opt_firstpass_fittingmethod;
    //
    if(m_opt_continuumcomponent=="tplfit"){
        linemodel.m_opt_tplfit_dustFit = Int32(m_opt_tplfit_dustfit=="yes");
        linemodel.m_opt_tplfit_extinction = Int32(m_opt_tplfit_igmfit=="yes");
        linemodel.m_opt_fitcontinuum_maxN = m_opt_continuumfitcount;
        linemodel.m_opt_tplfit_ignoreLinesSupport = Int32(m_opt_tplfit_ignoreLinesSupport=="yes");
        linemodel.m_opt_secondpasslcfittingmethod = m_opt_secondpasslcfittingmethod;

    }

    if(m_opt_rigidity=="tplshape")
    {
        linemodel.m_opt_tplratio_ismFit = Int32(m_opt_tplratio_ismfit=="yes");
        linemodel.m_opt_firstpass_tplratio_ismFit = Int32(m_opt_firstpass_tplratio_ismfit=="yes");
    }

    if(m_opt_continuumcomponent=="fromspectrum"){
        //check the continuum validity
        if( !_spcContinuum.IsFluxValid( lambdaRange.GetBegin(), lambdaRange.GetEnd() ) ){
            Log.LogWarning("Line Model - Failed to validate continuum spectrum flux on wavelength range (%.1f ; %.1f)",lambdaRange.GetBegin(), lambdaRange.GetEnd() );
            //throw std::runtime_error("Failed to validate continuum  flux");
        }else{
            Log.LogDetail( "Line Model - Successfully validated continuum flux on wavelength range (%.1f ; %.1f)", lambdaRange.GetBegin(), lambdaRange.GetEnd() );
        }

        /*
        //check the continuum flux axis
        const CSpectrumSpectralAxis& contspectralAxis = _spcContinuum.GetSpectralAxis();
        const CSpectrumFluxAxis& contfluxAxis = _spcContinuum.GetFluxAxis();
        for(Int32 j=0; j<contfluxAxis.GetSamplesCount(); j++)
        {
            if(isnan(contfluxAxis[j]))
            {
                Log.LogError( "Linemodelsolve(): NaN value found for the ContinuumFluxAxis at lambda=%f", contspectralAxis[j] );
                throw runtime_error("Linemodelsolve() NaN value found");
                break;
            }
        }
        //*/
    }

    //**************************************************
    //FIRST PASS
    //**************************************************
    Int32 retFirstPass = linemodel.ComputeFirstPass(dataStore,
                                          _spc,
                                          _spcContinuum,
                                          tplCatalog,
                                          tplCategoryList,
                                          m_calibrationPath,
                                          restraycatalog,
                                          m_opt_linetypefilter,
                                          m_opt_lineforcefilter,
                                          lambdaRange,
                                          m_opt_fittingmethod,
                                          m_opt_continuumcomponent,
                                          m_opt_lineWidthType,
                                          m_opt_resolution,
                                          m_opt_velocity_emission,
                                          m_opt_velocity_absorption,
                                          m_opt_continuumreest,
                                          m_opt_rules,
                                          m_opt_velocityfit,
                                          m_opt_firstpass_largegridstep,
                                          m_opt_firstpass_largegridsampling,
                                          m_opt_rigidity,
                                          m_opt_tplratio_reldirpath,
                                          m_opt_offsets_reldirpath);
    if( retFirstPass!=0 )
    {
        Log.LogError( "Line Model, first pass failed. Aborting" );
        return false;
    }

    //**************************************************
    //Compute z-candidates
    //**************************************************
    Bool overrideUseBestchi2forCandidates = false;
    // Int32 sign = 1;
    std::vector<Float64> fvals;
    std::shared_ptr<const CLineModelResult> lmresult = std::dynamic_pointer_cast<const CLineModelResult>( linemodel.getResult() );
    if(overrideUseBestchi2forCandidates)
    {
      // sign = -1;
        fvals = lmresult->ChiSquare;
    }else{
        std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
        std::shared_ptr<CPdfLogResult> zpriorResult = std::shared_ptr<CPdfLogResult>(new CPdfLogResult());
        Int32 retCombinePdf = CombinePDF(lmresult, m_opt_rigidity, m_opt_pdfcombination, m_opt_stronglinesprior, m_opt_euclidNHaEmittersPriorStrength, postmargZResult, zpriorResult);

        if(retCombinePdf!=0)
        {
            Log.LogError("Linemodel: Candidates search - Pdfz computation failed");
            return false;
        }else{
          // sign = -1;
            fvals = postmargZResult->valProbaLog;
        }
    }
    Int32 retCandidates = linemodel.ComputeCandidates(m_opt_extremacount, 1, fvals, m_opt_candidatesLogprobaCutThreshold); //BUG ? sign should change if chi2 or pdf are used... ?
    if( retCandidates!=0 )
    {
        Log.LogError( "Linemodel: Search for z-candidates failed. Aborting" );
        throw std::runtime_error("Linemodel: Search for z-candidates failed. Aborting");
    }


    //**************************************************
    //SECOND PASS
    //**************************************************
    bool skipSecondPass = (m_opt_skipsecondpass=="yes");
    if(!skipSecondPass)
    {
        Int32 retSecondPass = linemodel.ComputeSecondPass(dataStore,
                                                          _spc,
                                                          _spcContinuum,
                                                          tplCatalog,
                                                          tplCategoryList,
                                                          m_calibrationPath,
                                                          restraycatalog,
                                                          m_opt_linetypefilter,
                                                          m_opt_lineforcefilter,
                                                          lambdaRange,
                                                          m_opt_extremacount,
                                                          m_opt_fittingmethod,
                                                          m_opt_continuumcomponent,
                                                          m_opt_lineWidthType,
                                                          m_opt_resolution,
                                                          m_opt_velocity_emission,
                                                          m_opt_velocity_absorption,
                                                          m_opt_continuumreest,
                                                          m_opt_rules,
                                                          m_opt_velocityfit,
                                                          m_opt_rigidity,
                                                          m_opt_em_velocity_fit_min,
                                                          m_opt_em_velocity_fit_max,
                                                          m_opt_em_velocity_fit_step,
                                                          m_opt_abs_velocity_fit_min,
                                                          m_opt_abs_velocity_fit_max,
                                                          m_opt_abs_velocity_fit_step);
        if( retSecondPass!=0 )
        {
            Log.LogError( "Line Model, second pass failed. Aborting" );
            return false;
        }
    }
    std::shared_ptr<const CLineModelResult> result = std::dynamic_pointer_cast<const CLineModelResult>( linemodel.getResult() );


    if( !result )
    {
        //Log.LogInfo( "Failed to compute linemodel");
        return false;
    }else{
        //save linemodel chisquare results
        dataStore.StoreScopedGlobalResult( scopeStr.c_str(), result );
        //save linemodel extrema results
        std::string extremaResultsStr=scopeStr.c_str();
        extremaResultsStr.append("_extrema");
        //Log.LogError("Line Model, saving extrema results: %s", extremaResultsStr.c_str());
        dataStore.StoreScopedGlobalResult( extremaResultsStr.c_str(), result->GetExtremaResult() );

        //save linemodel firstpass extrema results
        std::string firstpassExtremaResultsStr=scopeStr.c_str();
        firstpassExtremaResultsStr.append("_firstpass_extrema");
        //Log.LogError("Line Model, saving firstpass extrema results: %s", firstpassExtremaResultsStr.c_str());
        dataStore.StoreScopedGlobalResult( firstpassExtremaResultsStr.c_str(), linemodel.GetFirstpassExtremaResult() );

        //save linemodel fitting and spectrum-model results
        linemodel.storeGlobalModelResults(dataStore);
    }

    return true;
}


#include <RedshiftLibrary/method/linemodelsolve.h>

#include <RedshiftLibrary/log/log.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/statistics/zprior.h>
#include <RedshiftLibrary/statistics/deltaz.h>
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
CLineModelSolve::CLineModelSolve(string calibrationPath):
    m_calibrationPath(calibrationPath)
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
    desc.append("\tparam: linemodel.continuumfit.method = {""templatefitting"", ""templatefittinglog""}\n");
    desc.append("\tparam: linemodel.continuumfit.ismfit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.continuumfit.igmfit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.continuumfit.count = <float value>\n");
    desc.append("\tparam: linemodel.continuumfit.ignorelinesupport = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.continuumfit.priors.betaA = <float value>\n");
    desc.append("\tparam: linemodel.continuumfit.priors.betaTE = <float value>\n");
    desc.append("\tparam: linemodel.continuumfit.priors.betaZ = <float value>\n");
    desc.append("\tparam: linemodel.continuumfit.priors.catalog_dirpath = <path>\n");
    desc.append("\tparam: linemodel.secondpasslcfittingmethod = {""no"", ""svdlcp2""}\n");
    desc.append("\tparam: linemodel.rigidity = {""rules"", ""tplcorr"", ""tplshape""}\n");
    desc.append("\tparam: linemodel.tplratio_catalog = <relative path>\n");
    desc.append("\tparam: linemodel.tplratio_ismfit = {""no"", ""yes""}\n");

    desc.append("\tparam: linemodel.tplratio.priors.betaA = <float value>\n");
    desc.append("\tparam: linemodel.tplratio.priors.betaTE = <float value>\n");
    desc.append("\tparam: linemodel.tplratio.priors.betaZ = <float value>\n");
    desc.append("\tparam: linemodel.tplratio.priors.catalog_dirpath = <path>\n");

    desc.append("\tparam: linemodel.offsets_catalog = <relative path>\n");
    desc.append("\tparam: linemodel.enableLSF = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.linewidthtype = {""instrumentdriven"", ""velocitydriven"",  ""combined"",  ""nispvsspsf201707"", ""fixed""}\n");
    desc.append("\tparam: linemodel.nsigmasupport = <float value>\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.rules = {""all"", ""balmer"", ""strongweak"", ""superstrong"", ""ratiorange"", ""ciiiratio"", ""no""}\n");
    desc.append("\tparam: linemodel.improveBalmerFit = {""yes"", ""no""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");
    desc.append("\tparam: linemodel.extremacountB = <float value>\n");
    desc.append("\tparam: linemodel.extremacutprobathreshold = <float value> (-1:disabled)\n");
    desc.append("\tparam: linemodel.velocityfit = {""yes"", ""no""}\n");
    desc.append("\tparam: linemodel.emvelocityfitmin = <float value>\n");
    desc.append("\tparam: linemodel.emvelocityfitmax = <float value>\n");
    desc.append("\tparam: linemodel.emvelocityfitstep = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitmin = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitmax = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitstep = <float value>\n");

    desc.append("\tparam: linemodel.lyaforcefit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.lyaforcedisablefit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.lyafit.asymfitmin = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.asymfitmax = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.asymfitstep = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.widthcoefffitmin = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.widthcoefffitmax = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.widthcoefffitstep = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.deltafitmin = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.deltafitmax = <float value>\n");
    desc.append("\tparam: linemodel.lyafit.deltafitstep = <float value>\n");

    //first pass
    desc.append("\tparam: linemodel.firstpass.largegridstep = <float value>, deactivated if negative or zero\n");
    desc.append("\tparam: linemodel.firstpass.tplratio_ismfit = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.firstpass.multiplecontinuumfit_disable = {""no"", ""yes""}\n");

    //second pass
    desc.append("\tparam: linemodel.skipsecondpass = {""no"", ""yes""}\n");
    desc.append("\tparam: linemodel.secondpass.continuumfit = {""fromfirstpass"", ""retryall"", ""refitfirstpass""}\n");
    desc.append("\tparam: linemodel.secondpass.halfwindowsize = <float value>\n");

    desc.append("\tparam: linemodel.pdfcombination = {""marg"", ""bestchi2""}\n");
    desc.append("\tparam: linemodel.pdf.margampcorr = {""yes"", ""no""}\n");

    desc.append("\tparam: linemodel.stronglinesprior = <float value>, penalization factor = positive value or -1 to deactivate\n");
    desc.append("\tparam: linemodel.haprior = <float value>, penalization factor = positive value (typical 1e-1 to 1e-5) or -1 to deactivate\n");
    desc.append("\tparam: linemodel.euclidnhaemittersStrength = <float value>, prior strength factor = positive value (typically 1 to 5) or -1 to deactivate\n");
    desc.append("\tparam: linemodel.modelpriorzStrength = <float value>, prior strength factor = positive value (typically 1 to 5) or -1 to deactivate\n");

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
    dataStore.GetScopedParam( "linemodel.secondpass.continuumfit", m_opt_secondpass_continuumfit, "fromfirstpass" );
    dataStore.GetScopedParam( "linemodel.secondpass.halfwindowsize", m_opt_secondpass_halfwindowsize, 0.005 );
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
    if(m_opt_continuumcomponent=="tplfit" || m_opt_continuumcomponent == "tplfitauto"){
        dataStore.GetScopedParam( "linemodel.continuumfit.method", m_opt_tplfit_method, "templatefittinglog" );
        m_opt_tplfit_method_secondpass = m_opt_tplfit_method; //for the moment, we use the same method for both passes
        dataStore.GetScopedParam( "linemodel.continuumfit.ismfit", m_opt_tplfit_dustfit, "yes" );
        dataStore.GetScopedParam( "linemodel.continuumfit.igmfit", m_opt_tplfit_igmfit, "yes" );
        dataStore.GetScopedParam( "linemodel.continuumfit.count", m_opt_continuumfitcount, 1 );
        dataStore.GetScopedParam( "linemodel.continuumfit.ignorelinesupport", m_opt_tplfit_ignoreLinesSupport, "no" );
        dataStore.GetScopedParam( "linemodel.continuumfit.priors.betaA", m_opt_tplfit_continuumprior_betaA, 1. );
        dataStore.GetScopedParam( "linemodel.continuumfit.priors.betaTE", m_opt_tplfit_continuumprior_betaTE, 1. );
        dataStore.GetScopedParam( "linemodel.continuumfit.priors.betaZ", m_opt_tplfit_continuumprior_betaZ, 1. );
        dataStore.GetScopedParam( "linemodel.continuumfit.priors.catalog_dirpath", m_opt_tplfit_continuumprior_dirpath, "" ); //no priors by default
    }
    dataStore.GetScopedParam( "linemodel.rigidity", m_opt_rigidity, "rules" );

    if(m_opt_rigidity=="tplshape")
    {
        dataStore.GetScopedParam( "linemodel.tplratio_catalog", m_opt_tplratio_reldirpath, "linecatalogs_tplshapes/linecatalogs_tplshape_ExtendedTemplatesJan2017v3_20170602_B14C_v5_emission" );
        dataStore.GetScopedParam( "linemodel.tplratio_ismfit", m_opt_tplratio_ismfit, "yes" );

        dataStore.GetScopedParam( "linemodel.tplratio.priors.betaA", m_opt_tplratio_prior_betaA, 1. );
        dataStore.GetScopedParam( "linemodel.tplratio.priors.betaTE", m_opt_tplratio_prior_betaTE, 1. );
        dataStore.GetScopedParam( "linemodel.tplratio.priors.betaZ", m_opt_tplratio_prior_betaZ, 1. );
        dataStore.GetScopedParam( "linemodel.tplratio.priors.catalog_dirpath", m_opt_tplratio_prior_dirpath, "" ); //no priors by default
    }else if(m_opt_rigidity=="rules")
    {
        dataStore.GetScopedParam( "linemodel.improveBalmerFit", m_opt_enableImproveBalmerFit, "yes" );
    }
    dataStore.GetScopedParam( "linemodel.offsets_catalog", m_opt_offsets_reldirpath, "linecatalogs_offsets/offsetsCatalogs_20170410_m150" );

    dataStore.GetScopedParam( "linemodel.enableLSF", m_opt_enableLSF, "no" );
    dataStore.GetScopedParam( "linemodel.linewidthtype", m_opt_lineWidthType, "velocitydriven" );
    dataStore.GetScopedParam( "linemodel.nsigmasupport", m_opt_nsigmasupport, 8.0 );
    dataStore.GetScopedParam( "linemodel.instrumentresolution", m_opt_resolution, 2350.0 );
    dataStore.GetScopedParam( "linemodel.velocityemission", m_opt_velocity_emission, 200.0 );
    dataStore.GetScopedParam( "linemodel.velocityabsorption", m_opt_velocity_absorption, 300.0 );
    dataStore.GetScopedParam( "linemodel.velocityfit", m_opt_velocityfit, "yes" );
    if(m_opt_velocityfit=="yes"){
        dataStore.GetScopedParam( "linemodel.emvelocityfitmin", m_opt_em_velocity_fit_min, 20.0 );
        dataStore.GetScopedParam( "linemodel.emvelocityfitmax", m_opt_em_velocity_fit_max, 300.0 );
        dataStore.GetScopedParam( "linemodel.emvelocityfitstep", m_opt_em_velocity_fit_step, 20.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitmin", m_opt_abs_velocity_fit_min, 150.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitmax", m_opt_abs_velocity_fit_max, 500.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitstep", m_opt_abs_velocity_fit_step, 50.0 );
    }
    dataStore.GetScopedParam( "linemodel.lyaforcefit", m_opt_lya_forcefit, "no" );
    dataStore.GetScopedParam( "linemodel.lyaforcedisablefit", m_opt_lya_forcedisablefit, "no" );
    dataStore.GetScopedParam( "linemodel.lyafit.asymfitmin", m_opt_lya_fit_asym_min, 0.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.asymfitmax", m_opt_lya_fit_asym_max, 4.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.asymfitstep", m_opt_lya_fit_asym_step, 1.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.widthfitmin", m_opt_lya_fit_width_min, 1.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.widthfitmax", m_opt_lya_fit_width_max, 4.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.widthfitstep", m_opt_lya_fit_width_step, 1.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.deltafitmin", m_opt_lya_fit_delta_min, 0.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.deltafitmax", m_opt_lya_fit_delta_max, 0.0 );
    dataStore.GetScopedParam( "linemodel.lyafit.deltafitstep", m_opt_lya_fit_delta_step, 1.0 );

    dataStore.GetScopedParam( "linemodel.continuumreestimation", m_opt_continuumreest, "no" );
    dataStore.GetScopedParam( "linemodel.rules", m_opt_rules, "all" );
    dataStore.GetScopedParam( "linemodel.extremacount", m_opt_extremacount, 10 );
    dataStore.GetScopedParam( "linemodel.extremacountB", m_opt_extremacountB, 0 );
    dataStore.GetScopedParam( "linemodel.extremacutprobathreshold", m_opt_candidatesLogprobaCutThreshold, -1 );
    dataStore.GetScopedParam( "linemodel.stronglinesprior", m_opt_stronglinesprior, -1);
    dataStore.GetScopedParam( "linemodel.haprior", m_opt_haPrior, -1);
    dataStore.GetScopedParam( "linemodel.euclidnhaemittersStrength", m_opt_euclidNHaEmittersPriorStrength, -1);
    dataStore.GetScopedParam( "linemodel.modelpriorzStrength", m_opt_modelZPriorStrength, -1);
    dataStore.GetScopedParam( "linemodel.pdfcombination", m_opt_pdfcombination, "marg");
    dataStore.GetScopedParam( "linemodel.pdf.margampcorr", m_opt_pdf_margAmpCorrection, "no");
    dataStore.GetScopedParam( "linemodel.saveintermediateresults", m_opt_saveintermediateresults, "no");

    //Auto-correct fitting method
    std::string forcefittingmethod = "individual";
    if(m_opt_rigidity=="tplshape" && m_opt_fittingmethod == "hybrid")
    {
        m_opt_fittingmethod = forcefittingmethod;
        dataStore.SetScopedParam("linemodel.fittingmethod", m_opt_fittingmethod);
        Log.LogInfo( "Linemodel fitting method auto-correct due to tplshape rigidity");
    }
    if(m_opt_rigidity=="tplshape" && m_opt_firstpass_fittingmethod == "hybrid")
    {
        m_opt_firstpass_fittingmethod = forcefittingmethod;
        dataStore.SetScopedParam("linemodel.firstpass.fittingmethod", m_opt_firstpass_fittingmethod);
        Log.LogInfo( "Linemodel first pass fitting method auto-correct due to tplshape rigidity");
    }

    Log.LogInfo( "Linemodel parameters:");
    Log.LogInfo( "    -linetypefilter: %s", m_opt_linetypefilter.c_str());
    Log.LogInfo( "    -lineforcefilter: %s", m_opt_lineforcefilter.c_str());
    Log.LogInfo( "    -fittingmethod: %s", m_opt_fittingmethod.c_str());
    Log.LogInfo( "    -enableLSF: %s", m_opt_enableLSF.c_str());
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

    Log.LogInfo( "    -nsigmasupport: %.1f", m_opt_nsigmasupport);

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
        if(m_opt_fittingmethod == "hybrid")
        {
            Log.LogInfo( "      -Improve Balmer Fit: %s", m_opt_enableImproveBalmerFit.c_str());
        }
    }else if(m_opt_rigidity=="tplshape")
    {
        Log.LogInfo( "      -tplratio_catalog: %s", m_opt_tplratio_reldirpath.c_str());
        Log.LogInfo( "      -tplratio_ismfit: %s", m_opt_tplratio_ismfit.c_str());
        Log.LogInfo( "      -tplfit_priors_dirpath: %s", m_opt_tplratio_prior_dirpath.c_str());
        Log.LogInfo( "      -tplratio_priors_betaA:  %f", m_opt_tplratio_prior_betaA);
        Log.LogInfo( "      -tplratio_priors_betaTE:  %f", m_opt_tplratio_prior_betaTE);
        Log.LogInfo( "      -tplratio_priors_betaZ:  %f", m_opt_tplratio_prior_betaZ);
    }
    Log.LogInfo( "    -linemodel offsets_catalog: %s", m_opt_offsets_reldirpath.c_str());

    Log.LogInfo( "    -continuumcomponent: %s", m_opt_continuumcomponent.c_str());
    if(m_opt_continuumcomponent=="tplfit" || m_opt_continuumcomponent=="tplfitauto"){
        Log.LogInfo( "      -tplfit_method: %s", m_opt_tplfit_method.c_str());
        Log.LogInfo( "      -tplfit_ismfit: %s", m_opt_tplfit_dustfit.c_str());
        Log.LogInfo( "      -tplfit_igmfit: %s", m_opt_tplfit_igmfit.c_str());
        Log.LogInfo( "      -continuum fit count:  %.0f", m_opt_continuumfitcount);
        Log.LogInfo( "      -tplfit_ignorelinesupport: %s", m_opt_tplfit_ignoreLinesSupport.c_str());
        Log.LogInfo( "      -tplfit_secondpass-LC-fitting-method: %s", m_opt_secondpasslcfittingmethod.c_str());
        Log.LogInfo( "      -tplfit_priors_dirpath: %s", m_opt_tplfit_continuumprior_dirpath.c_str());
        Log.LogInfo( "      -tplfit_priors_betaA:  %f", m_opt_tplfit_continuumprior_betaA);
        Log.LogInfo( "      -tplfit_priors_betaTE:  %f", m_opt_tplfit_continuumprior_betaTE);
        Log.LogInfo( "      -tplfit_priors_betaZ:  %f", m_opt_tplfit_continuumprior_betaZ);
    }
    Log.LogInfo( "    -continuumreestimation: %s", m_opt_continuumreest.c_str());
    Log.LogInfo( "    -extremacount: %.0f", m_opt_extremacount);
    Log.LogInfo( "    -extremacount-firstpass B: %.0f", m_opt_extremacountB);
    Log.LogInfo( "    -extrema cut proba-threshold: %.0f", m_opt_candidatesLogprobaCutThreshold);
    Log.LogInfo( "    -first pass:");
    Log.LogInfo( "      -largegridstep: %.6f", m_opt_firstpass_largegridstep);
    Log.LogInfo( "      -fittingmethod: %s", m_opt_firstpass_fittingmethod.c_str());
    Log.LogInfo( "      -tplratio_ismfit: %s", m_opt_firstpass_tplratio_ismfit.c_str());
    Log.LogInfo( "      -multiplecontinuumfit_disable: %s", m_opt_firstpass_disablemultiplecontinuumfit.c_str());

    Log.LogInfo( "    -second pass:");
    Log.LogInfo( "      -skip second pass: %s", m_opt_skipsecondpass.c_str());
    Log.LogInfo( "      -continuum fit method: %s", m_opt_secondpass_continuumfit.c_str());

    Log.LogInfo( "    -pdf-stronglinesprior: %e", m_opt_stronglinesprior);
    Log.LogInfo( "    -pdf-hapriorstrength: %e", m_opt_haPrior);
    Log.LogInfo( "    -pdf-euclidNHaEmittersPriorStrength: %e", m_opt_euclidNHaEmittersPriorStrength);
    Log.LogInfo( "    -pdf-modelpriorzStrength: %e", m_opt_modelZPriorStrength);
    Log.LogInfo( "    -pdf-combination: %s", m_opt_pdfcombination.c_str()); // "marg";    // "bestchi2";    // "bestproba";
    Log.LogInfo( "    -pdf-margAmpCorrection: %s", m_opt_pdf_margAmpCorrection.c_str());

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
                                                                 const CTemplateCatalog& tplCatalog,
                                                                 const TStringList& tplCategoryList,
                                                                 const CRayCatalog& restraycatalog,
                                                                 const TFloat64Range& lambdaRange,
                                                                 const TFloat64List& redshifts,
                                                                 const std::string outputPdfRelDir,
                                                                 const Float64 radius )
{
    CDataStore::CAutoScope resultScope( dataStore, "linemodelsolve" );
    m_outputPdfRelDir = outputPdfRelDir;
    m_redshiftSeparation = radius;

    PopulateParameters( dataStore );


    /* ------------------------  Solve the Linemodel  --------------------------  */

    Int32 retSolve = Solve( dataStore, spc, tplCatalog, tplCategoryList, restraycatalog, lambdaRange, redshifts );

    if(!retSolve){
        return NULL;
    }

    std::string scope = dataStore.GetCurrentScopeName() + ".linemodel";
    auto results = dataStore.GetGlobalResult( scope.c_str() );
    if(results.expired())
    {
        Log.LogError("linemodelsolve: Unable to retrieve linemodel results");
        return NULL;
    }
    std::shared_ptr<const CLineModelResult> result = std::dynamic_pointer_cast<const CLineModelResult>( results.lock());

    //Save chisquareTplshape results
    if(m_opt_enableSaveChisquareTplshapeResults)
        StoreChisquareTplShapeResults(dataStore, result);

    //std::shared_ptr<CPdfLogResult> zpriorResult = std::make_shared<CPdfLogResult>();

    // prepare the linemodel chisquares and prior results for pdf computation
    ChisquareArray chisquares = BuildChisquareArray(result,
                                                    m_opt_rigidity,
                                                    m_opt_pdfcombination,
                                                    m_opt_stronglinesprior,
                                                    m_opt_haPrior,
                                                    m_opt_euclidNHaEmittersPriorStrength,
                                                    m_opt_modelZPriorStrength);
                                                    
    /*
    Log.LogDetail("    linemodelsolve: Storing priors (size=%d)", zpriorResult->Redshifts.size());
    std::string priorPath = outputPdfRelDir+"/logprior.logP_Z_data";
    dataStore.StoreGlobalResult( priorPath.c_str(), zpriorResult);
    */

    /* ------------------------  COMPUTE POSTERIOR PDF  --------------------------  */

    COperatorPdfz pdfz(m_opt_pdfcombination, 
                        0.0, // no peak Separation in 2nd pass
                        0.0, // cut threshold
                        m_opt_extremacount, // max nb of final (2nd pass) candidates
                        "SPE", //Id_prefix
                        false, // do not allow extrema at border
                        1,  // one peak/window only
                        m_linemodel.m_secondpass_parameters_extremaResult.ExtendedRedshifts,
                        m_linemodel.m_secondpass_parameters_extremaResult.GetIDs()
                        ); 
    
    std::shared_ptr<CPdfCandidateszResult> candidateResult = pdfz.Compute(chisquares);

    // store PDF results
    Log.LogInfo("%s: Storing PDF results", __func__);
    std::string pdfPath = outputPdfRelDir+"/logposterior.logMargP_Z_data";
    dataStore.StoreGlobalResult( pdfPath.c_str(), pdfz.m_postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
    dataStore.StoreGlobalResult("candidatesresult", candidateResult);

    // Get linemodel results at extrema (recompute spectrum model etc.)
    std::shared_ptr<const CLineModelExtremaResult> ExtremaResult =
        m_linemodel.SaveExtremaResults( spc, lambdaRange, candidateResult->m_ranked_candidates, 
                                        m_opt_continuumreest);

    // store extrema results
    storeExtremaResults(dataStore, ExtremaResult );

    //SaveContinuumPDF(dataStore, result);
    // TBD

    // create the solveresult
    std::shared_ptr<CLineModelSolveResult> lmsolveresult = 
        std::make_shared<CLineModelSolveResult>( ExtremaResult, 
                                                 m_opt_pdfcombination,
                                                 pdfz.m_postmargZResult->valEvidenceLog);


    return lmsolveresult;
}

ChisquareArray CLineModelSolve::BuildChisquareArray(std::shared_ptr<const CLineModelResult> result,
                                                    std::string opt_rigidity,
                                                    std::string opt_combine,
                                                    Float64 opt_stronglinesprior,
                                                    Float64 opt_hapriorstrength,
                                                    Float64 opt_euclidNHaEmittersPriorStrength,
                                                    Float64 opt_modelPriorZStrength) const
{
    Log.LogDetail("LinemodelSolve: building chisquare array");

    ChisquareArray chisquarearray;
    
    chisquarearray.cstLog = result->cstLog;
    Log.LogDetail("%s: using cstLog = %f", __func__, chisquarearray.cstLog);
        
    chisquarearray.redshifts = result->Redshifts; 

    bool zPriorStrongLinePresence = (opt_stronglinesprior>0.0);
    if(zPriorStrongLinePresence)
    {
        Log.LogDetail("%s: StrongLinePresence prior enabled: factor=%e", __func__, opt_stronglinesprior);
    }else{
        Log.LogDetail("%s: StrongLinePresence prior disabled", __func__);
    }
    bool zPriorHaStrongestLine = (opt_hapriorstrength>0.0);
    if(zPriorHaStrongestLine)
    {
        Log.LogDetail("%s: Ha strongest line prior enabled: factor=%e", __func__, opt_hapriorstrength);
    }else{
        Log.LogDetail("%s: Ha strongest line prior disabled", __func__);
    }

    Float64 opt_nlines_snr_penalization_factor = -1;
    bool zPriorNLineSNR = (opt_nlines_snr_penalization_factor>0.0);
    if(zPriorNLineSNR)
    {
        Log.LogDetail("%s: N lines snr>cut prior enabled: factor=%e", __func__, opt_nlines_snr_penalization_factor);
    }else{
        Log.LogDetail("%s: N lines snr>cut prior disabled", __func__);
    }

    //hardcoded Euclid-NHaZprior parameter
    bool zPriorEuclidNHa = false;
    if(opt_euclidNHaEmittersPriorStrength>0.0)
    {
        zPriorEuclidNHa = true;
        Log.LogDetail("%s: EuclidNHa prior enabled, with strength-coeff: %e", __func__, opt_euclidNHaEmittersPriorStrength);
    }else{
        Log.LogDetail("%s: EuclidNHa prior disabled", __func__);
    }

    bool zPriorLines = true;
    Log.LogDetail("%s: PriorLinesTplshapes.size()=%d", __func__, result->PriorLinesTplshapes.size());
    if( !boost::filesystem::exists( m_opt_tplratio_prior_dirpath ) || result->PriorLinesTplshapes.size()!=result->ChiSquareTplshapes.size())
    {
        zPriorLines = false;
    }
    if(zPriorLines)
    {
        Log.LogDetail("%s: Lines Prior enabled", __func__);
    }else{
        Log.LogDetail("%s: Lines Prior disabled", __func__);
    }

    CZPrior zpriorhelper;

    if(opt_rigidity!="tplshape" || opt_combine=="bestchi2")
    {
        chisquarearray.zpriors.emplace_back();
        TFloat64List & zpriors = chisquarearray.zpriors.back();

        if(zPriorStrongLinePresence)
        {
            UInt32 lineTypeFilter = 1;// for emission lines only
            TBoolList strongLinePresence = result->GetStrongLinesPresence(lineTypeFilter, result->LineModelSolutions);

            zpriors = zpriorhelper.GetStrongLinePresenceLogZPrior(strongLinePresence, opt_stronglinesprior);
        }else{
            zpriors = zpriorhelper.GetConstantLogZPrior(result->Redshifts.size());
        }
        if(zPriorHaStrongestLine)
        {
            TBoolList wHaStronglinePresence = result->GetStrongestLineIsHa(result->LineModelSolutions); //whasp for lm-tplratio
            std::vector<Float64> zlogPriorHaStrongest = zpriorhelper.GetStrongLinePresenceLogZPrior(wHaStronglinePresence, opt_hapriorstrength);
            zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorHaStrongest);
        }
        if(zPriorEuclidNHa)
        {
            std::vector<Float64> zlogPriorNHa = zpriorhelper.GetEuclidNhaLogZPrior(result->Redshifts, opt_euclidNHaEmittersPriorStrength);
            zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorNHa);
        }
        if(zPriorNLineSNR)
        {
            std::vector<Int32> n_lines_above_snr = result->GetNLinesAboveSnrcut(result->LineModelSolutions);
            std::vector<Float64> zlogPriorNLinesAboveSNR = zpriorhelper.GetNLinesSNRAboveCutLogZPrior(n_lines_above_snr, opt_nlines_snr_penalization_factor);
            zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorNLinesAboveSNR);
        }

        //correct chi2 if necessary
        TFloat64List logLikelihoodCorrected(result->ChiSquare.size(), DBL_MAX);
        for ( UInt32 k=0; k<result->Redshifts.size(); k++ )
        {
            logLikelihoodCorrected[k] = result->ChiSquare[k];
        }
        if(false && m_opt_pdf_margAmpCorrection=="yes") //maybe there should not be a scalemarg correction for the bestchi2 option ? Todo: raise warning then...
        {
            for ( UInt32 k=0; k<result->Redshifts.size(); k++ )
            {
                logLikelihoodCorrected[k] += result->ScaleMargCorrection[k];
            }
        }
        
        chisquarearray.chisquares.push_back(std::move(logLikelihoodCorrected));

    }else if(opt_combine=="bestproba" || opt_combine=="marg"){

        //Log.LogInfo("Linemodel: Pdfz computation - combination: method=%s, n=%d", opt_combine.c_str(), result->ChiSquareTplshapes.size());
        for(Int32 k=0; k<result->ChiSquareTplshapes.size(); k++)
        {
            chisquarearray.zpriors.emplace_back();
            TFloat64List & zpriors = chisquarearray.zpriors.back();
    
            if(zPriorStrongLinePresence)
            {
                TBoolList const & strongLinePresence = result->StrongELPresentTplshapes[k];
                zpriors = zpriorhelper.GetStrongLinePresenceLogZPrior(strongLinePresence, opt_stronglinesprior);
            }else
            {
                zpriors = zpriorhelper.GetConstantLogZPrior(result->Redshifts.size());
            }

            if(zPriorHaStrongestLine)
            {
                TBoolList wHaStronglinePresence = result->GetStrongestLineIsHa(result->LineModelSolutions); //whasp for lm-tplratio
                std::vector<Float64> zlogPriorHaStrongest = zpriorhelper.GetStrongLinePresenceLogZPrior(wHaStronglinePresence, opt_hapriorstrength);
                zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorHaStrongest);
            }
            if(zPriorEuclidNHa)
            {
                std::vector<Float64> zlogPriorNHa = zpriorhelper.GetEuclidNhaLogZPrior(result->Redshifts, opt_euclidNHaEmittersPriorStrength);
                zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorNHa);
            }
            if(zPriorNLineSNR)
            {
                std::vector<Int32> n_lines_above_snr = result->NLinesAboveSNRTplshapes[k];
                std::vector<Float64> zlogPriorNLinesAboveSNR = zpriorhelper.GetNLinesSNRAboveCutLogZPrior(n_lines_above_snr, opt_nlines_snr_penalization_factor);
                zpriors = zpriorhelper.CombineLogZPrior(zpriors, zlogPriorNLinesAboveSNR);
            }
        }

        //correct chi2 if necessary
        std::vector<TFloat64List> ChiSquareTplshapesCorrected;
        for(Int32 k=0; k<result->ChiSquareTplshapes.size(); k++)
        {
            chisquarearray.chisquares.emplace_back(result->ChiSquareTplshapes[k].size(), DBL_MAX);
            TFloat64List & logLikelihoodCorrected = chisquarearray.chisquares.back();
            
            logLikelihoodCorrected = result->ChiSquareTplshapes[k];

            if(m_opt_pdf_margAmpCorrection=="yes") //nb: this is experimental.
            {
                //find max scalemargcorr
                //*
                Float64 maxscalemargcorr=-DBL_MAX;
                for ( UInt32 kz=0; kz<result->Redshifts.size(); kz++ )
                {
                    if(maxscalemargcorr < result->ScaleMargCorrectionTplshapes[k][kz])
                    {
                        maxscalemargcorr = result->ScaleMargCorrectionTplshapes[k][kz];
                    }
                }
                Log.LogDetail("%s: maxscalemargcorr= %e", __func__,  maxscalemargcorr);
                //*/
                for ( UInt32 kz=0; kz<result->Redshifts.size(); kz++ )
                {
                    if(result->ScaleMargCorrectionTplshapes[k][kz]!=0) //warning, this is experimental.
                    {
                        logLikelihoodCorrected[kz] += result->ScaleMargCorrectionTplshapes[k][kz] - maxscalemargcorr;
                    }
                }
            }
            if(zPriorLines && result->PriorLinesTplshapes[k].size()==result->Redshifts.size())
            {
                for ( UInt32 kz=0; kz<result->Redshifts.size(); kz++ )
                {
                    logLikelihoodCorrected[kz] += result->PriorLinesTplshapes[k][kz];
                }
            }
        }

        chisquarearray.modelpriors = result->PriorTplshapes;


        // todo : store priors for each tplshape model ?
        //  ->  in ::compute with ChisquareArray

    }else{
        Log.LogError("Linemodel: Unable to parse pdf combination method option");
    }

    return chisquarearray;

}

/*
Int32 CLineModelSolve::SaveContinuumPDF(CDataStore& store, std::shared_ptr<const CLineModelResult> result)
{
    Log.LogInfo("Linemodel: continuum Pdfz computation");
    std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
    COperatorPdfz pdfz;
    CZPrior zpriorhelper;
    Float64 cstLog = result->cstLog;
    TFloat64List logProba;
    Float64 logEvidence;

    Int32 retPdfz=-1;

    std::shared_ptr<CPdfLogResult> zPrior = std::shared_ptr<CPdfLogResult>(new CPdfLogResult());
    zPrior->SetSize(result->Redshifts.size());
    for ( UInt32 k=0; k<result->Redshifts.size(); k++ )
    {
        zPrior->Redshifts[k] = result->Redshifts[k];
    }

    zPrior->valProbaLog = zpriorhelper.GetConstantLogZPrior(result->Redshifts.size());


    //correct chi2 if necessary: todo add switch
    TFloat64List logLikelihoodCorrected(result->ChiSquareContinuum.size(), DBL_MAX);
    for ( UInt32 k=0; k<result->Redshifts.size(); k++ )
    {
        logLikelihoodCorrected[k] = result->ChiSquareContinuum[k];// + result->ScaleMargCorrectionContinuum[k];
    }
    retPdfz = pdfz.ComputePdf(logLikelihoodCorrected, result->Redshifts, cstLog, zPrior->valProbaLog, logProba, logEvidence);
    if(retPdfz==0){
        store.StoreGlobalResult( "zPDF/logpriorcontinuum.logP_Z_data", zPrior);

        postmargZResult->countTPL = result->Redshifts.size(); // assumed 1 model per z
        postmargZResult->Redshifts.resize(result->Redshifts.size());
        postmargZResult->valProbaLog.resize(result->Redshifts.size());
        for ( UInt32 k=0; k<result->Redshifts.size(); k++ )
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
*/


///
/// \brief COperatorLineModel::storeGlobalModelResults
/// stores the linemodel results as global results in the datastore
///
void CLineModelSolve::storeExtremaResults( CDataStore &dataStore,
                                           std::shared_ptr<const CLineModelExtremaResult> ExtremaResult) const 
{
    std::string extremaResultsStr = "linemodel_extrema";
    Log.LogInfo("Linemodel, saving extrema results: %s", extremaResultsStr.c_str());
    dataStore.StoreScopedGlobalResult( extremaResultsStr.c_str(), ExtremaResult );

    Int32 nResults = ExtremaResult->size();

    for (Int32 k = 0; k < nResults; k++)
    {
        std::string fname_spc =
            (boost::format("linemodel_spc_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(fname_spc.c_str(),
                                          ExtremaResult->m_savedModelSpectrumResults[k]);

        std::string fname_fit =
            (boost::format("linemodel_fit_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(fname_fit.c_str(),
                                          ExtremaResult->m_savedModelFittingResults[k]);

        std::string fname_fitcontinuum =
            (boost::format("linemodel_fitcontinuum_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(fname_fitcontinuum.c_str(), 
                                          ExtremaResult->m_savedModelContinuumFittingResults[k]);

        std::string fname_rules =
            (boost::format("linemodel_rules_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(fname_rules.c_str(),
                                          ExtremaResult->m_savedModelRulesResults[k]);

        std::string nameBaselineStr =
            (boost::format("linemodel_continuum_extrema_%1%") % k).str();
        dataStore.StoreScopedGlobalResult(nameBaselineStr.c_str(), 
                                         ExtremaResult->m_savedModelContinuumSpectrumResults[k]);
    }
}


void CLineModelSolve::StoreChisquareTplShapeResults(CDataStore & dataStore, std::shared_ptr<const CLineModelResult> result) const
{
    for(Int32 km=0; km<result->ChiSquareTplshapes.size(); km++)
    {
        std::shared_ptr<CLineModelResult> result_chisquaretplshape = std::shared_ptr<CLineModelResult>( new CLineModelResult() );
        result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0, std::vector<Float64>() );
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
        result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0, std::vector<Float64>() );
        for(Int32 kz=0; kz<result->Redshifts.size(); kz++)
        {
            result_chisquaretplshape->ChiSquare[kz] = result->ScaleMargCorrectionTplshapes[km][kz];
        }

        std::string resname = (boost::format("linemodel_chisquaretplshape/linemodel_scalemargcorrtplshape_%d") % km).str();
        dataStore.StoreScopedGlobalResult( resname.c_str(), result_chisquaretplshape );
    }

    //Save PriorLinesTplshapes results
    for(Int32 km=0; km<result->PriorLinesTplshapes.size(); km++)
    {
        std::shared_ptr<CLineModelResult> result_chisquaretplshape = std::shared_ptr<CLineModelResult>( new CLineModelResult() );
        result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0, std::vector<Float64>() );
        for(Int32 kz=0; kz<result->Redshifts.size(); kz++)
        {
            result_chisquaretplshape->ChiSquare[kz] = result->PriorLinesTplshapes[km][kz];
        }

        std::string resname = (boost::format("linemodel_chisquaretplshape/linemodel_priorlinestplshape_%d") % km).str();
        dataStore.StoreScopedGlobalResult( resname.c_str(), result_chisquaretplshape );
    }

    //Save PriorContinuumTplshapes results
    std::shared_ptr<CLineModelResult> result_chisquaretplshape = std::shared_ptr<CLineModelResult>( new CLineModelResult() );
    result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0, std::vector<Float64>() );
    for(Int32 kz=0; kz<result->Redshifts.size(); kz++)
    {
        result_chisquaretplshape->ChiSquare[kz] = result->ContinuumModelSolutions[kz].tplLogPrior;
    }

    std::string resname = (boost::format("linemodel_chisquaretplshape/linemodel_priorcontinuumtplshape")).str();
    dataStore.StoreScopedGlobalResult( resname.c_str(), result_chisquaretplshape );

}

/**
 * \brief
 * Retrieve the true-velocities from a hardcoded ref file path
 * nb: this is a hack for development purposes
 **/
Int32 getVelocitiesFromRefFile(const char* filePath, std::string spcid, Float64& elv, Float64& alv)
{
    std::ifstream file;

    file.open( filePath, std::ifstream::in );
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
 * Create a continuum object by subtracting spcWithoutContinuum from the spc.
 * Configure the opt_XXX variables from the dataStore scope parameters.
 * LogInfo the opt_XXX values.
 * Create a COperatorLineModel, call its Compute method.
 * If that returned true, store results.
 **/
Bool CLineModelSolve::Solve( CDataStore& dataStore,
                             const CSpectrum& spc,
                             const CTemplateCatalog& tplCatalog,
                             const TStringList& tplCategoryList,
                             const CRayCatalog& restraycatalog,
                             const TFloat64Range& lambdaRange,
                             const TFloat64List& redshifts )
{
    std::string scopeStr = "linemodel";

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
    Int32 retInit = m_linemodel.Init(spc, redshifts, m_opt_continuumcomponent, m_opt_nsigmasupport, m_opt_secondpass_halfwindowsize, m_redshiftSeparation);
    if( retInit!=0 )
    {
        Log.LogError( "Linemodel, init failed. Aborting" );
        throw std::runtime_error( "Linemodel, init failed. Aborting" );
    }
    m_linemodel.m_opt_firstpass_fittingmethod=m_opt_firstpass_fittingmethod;
    //
    if(m_opt_continuumcomponent=="tplfit" || m_opt_continuumcomponent=="tplfitauto"){
        m_linemodel.m_opt_tplfit_method = m_opt_tplfit_method;
        m_linemodel.m_opt_tplfit_method_secondpass = m_opt_tplfit_method_secondpass;
        m_linemodel.m_opt_tplfit_dustFit = Int32(m_opt_tplfit_dustfit=="yes");
        m_linemodel.m_opt_tplfit_extinction = Int32(m_opt_tplfit_igmfit=="yes");
        m_linemodel.m_opt_fitcontinuum_maxN = m_opt_continuumfitcount;
        m_linemodel.m_opt_tplfit_ignoreLinesSupport = Int32(m_opt_tplfit_ignoreLinesSupport=="yes");
        m_linemodel.m_opt_secondpasslcfittingmethod = m_opt_secondpasslcfittingmethod;
        m_linemodel.m_opt_tplfit_continuumprior_dirpath = m_opt_tplfit_continuumprior_dirpath;
        m_linemodel.m_opt_tplfit_continuumprior_betaA = m_opt_tplfit_continuumprior_betaA;
        m_linemodel.m_opt_tplfit_continuumprior_betaTE = m_opt_tplfit_continuumprior_betaTE;
        m_linemodel.m_opt_tplfit_continuumprior_betaZ = m_opt_tplfit_continuumprior_betaZ;
    }

    m_linemodel.m_opt_enableLSF=m_opt_enableLSF;

    m_linemodel.m_opt_lya_forcefit=m_opt_lya_forcefit;
    m_linemodel.m_opt_lya_forcedisablefit=m_opt_lya_forcedisablefit;
    m_linemodel.m_opt_lya_fit_asym_min=m_opt_lya_fit_asym_min;
    m_linemodel.m_opt_lya_fit_asym_max=m_opt_lya_fit_asym_max;
    m_linemodel.m_opt_lya_fit_asym_step=m_opt_lya_fit_asym_step;
    m_linemodel.m_opt_lya_fit_width_min=m_opt_lya_fit_width_min;
    m_linemodel.m_opt_lya_fit_width_max=m_opt_lya_fit_width_max;
    m_linemodel.m_opt_lya_fit_width_step=m_opt_lya_fit_width_step;
    m_linemodel.m_opt_lya_fit_delta_min=m_opt_lya_fit_delta_min;
    m_linemodel.m_opt_lya_fit_delta_max=m_opt_lya_fit_delta_max;
    m_linemodel.m_opt_lya_fit_delta_step=m_opt_lya_fit_delta_step;

    if(m_opt_rigidity=="tplshape")
    {
        m_linemodel.m_opt_tplratio_ismFit = Int32(m_opt_tplratio_ismfit=="yes");
        m_linemodel.m_opt_firstpass_tplratio_ismFit = Int32(m_opt_firstpass_tplratio_ismfit=="yes");

        m_linemodel.m_opt_tplratio_prior_dirpath = m_opt_tplratio_prior_dirpath;
        m_linemodel.m_opt_tplratio_prior_betaA = m_opt_tplratio_prior_betaA;
        m_linemodel.m_opt_tplratio_prior_betaTE = m_opt_tplratio_prior_betaTE;
        m_linemodel.m_opt_tplratio_prior_betaZ = m_opt_tplratio_prior_betaZ;
    }

    if(m_opt_rigidity=="rules")
    {
        m_linemodel.m_opt_enableImproveBalmerFit = m_opt_enableImproveBalmerFit;
    }

    //**************************************************
    //FIRST PASS
    //**************************************************
    Int32 retFirstPass = m_linemodel.ComputeFirstPass(spc,
                                                    tplCatalog,
                                                    tplCategoryList,
                                                    m_calibrationPath,
                                                    restraycatalog,
                                                    m_opt_linetypefilter,
                                                    m_opt_lineforcefilter,
                                                    lambdaRange,
                                                    m_opt_fittingmethod,
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
        Log.LogError( "Linemodel, first pass failed. Aborting" );
        throw runtime_error("Linemodel, first pass failed. Aborting");
        return false;
    }

    //**************************************************
    //Compute z-candidates
    //**************************************************
    std::shared_ptr<const CLineModelResult> lmresult = std::dynamic_pointer_cast<const CLineModelResult>( m_linemodel.getResult() );

    ChisquareArray chisquares = BuildChisquareArray(lmresult,
                                                    m_opt_rigidity,
                                                    m_opt_pdfcombination,
                                                    m_opt_stronglinesprior,
                                                    m_opt_haPrior,
                                                    m_opt_euclidNHaEmittersPriorStrength,
                                                    m_opt_modelZPriorStrength);

    //TODO deal with the case lmresult->Redshifts=1
    Int32 extremacount = 5;
    COperatorPdfz pdfz(m_opt_pdfcombination,
                        2*m_opt_secondpass_halfwindowsize, // peak separation
                        m_opt_candidatesLogprobaCutThreshold,
                        extremacount,
                        "FPE");

    std::shared_ptr<CPdfCandidateszResult> candResult = pdfz.Compute(chisquares, false);
    
    m_linemodel.SetFirstPassCandidates(candResult->m_ranked_candidates);
 
    //**************************************************
    //FIRST PASS + CANDIDATES - B
    //**************************************************
    Bool enableFirstpass_B = (m_opt_extremacountB>0) && (m_opt_continuumcomponent=="tplfit" || m_opt_continuumcomponent=="tplfitauto") && (m_opt_extremacountB>1);
    COperatorLineModel linemodel_fpb;
    std::string fpb_opt_continuumcomponent = "fromspectrum";//Note: this is hardocoded! given that condition for FPB relies on having "tplfit"
    Int32 retInitB = linemodel_fpb.Init(spc, redshifts, fpb_opt_continuumcomponent, m_opt_nsigmasupport, m_opt_secondpass_halfwindowsize, m_redshiftSeparation);
    if( retInitB!=0 )
    {
        Log.LogError( "Linemodel fpB, init failed. Aborting" );
        return false;
    }
    if(enableFirstpass_B)
    {
        Log.LogInfo( "Linemodel FIRST PASS B enabled. Computing now." );

        linemodel_fpb.m_opt_firstpass_fittingmethod=m_opt_firstpass_fittingmethod;
        //

        if(fpb_opt_continuumcomponent=="tplfit" || fpb_opt_continuumcomponent=="tplfitauto"){
            linemodel_fpb.m_opt_tplfit_dustFit = Int32(m_opt_tplfit_dustfit=="yes");
            linemodel_fpb.m_opt_tplfit_extinction = Int32(m_opt_tplfit_igmfit=="yes");
            linemodel_fpb.m_opt_fitcontinuum_maxN = m_opt_continuumfitcount;
            linemodel_fpb.m_opt_tplfit_ignoreLinesSupport = Int32(m_opt_tplfit_ignoreLinesSupport=="yes");
            linemodel_fpb.m_opt_secondpasslcfittingmethod = m_opt_secondpasslcfittingmethod;

        }

        if(m_opt_rigidity=="tplshape")
        {
            linemodel_fpb.m_opt_tplratio_ismFit = Int32(m_opt_tplratio_ismfit=="yes");
            linemodel_fpb.m_opt_firstpass_tplratio_ismFit = Int32(m_opt_firstpass_tplratio_ismfit=="yes");
        }

        //**************************************************
        //FIRST PASS B
        //**************************************************
        Int32 retFirstPass = linemodel_fpb.ComputeFirstPass(spc,
                                                            tplCatalog,
                                                            tplCategoryList,
                                                            m_calibrationPath,
                                                            restraycatalog,
                                                            m_opt_linetypefilter,
                                                            m_opt_lineforcefilter,
                                                            lambdaRange,
                                                            m_opt_fittingmethod,
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
            Log.LogError( "Linemodel, first pass failed. Aborting" );
            return false;
        }

        //**************************************************
        //Compute z-candidates B
        //**************************************************
        std::shared_ptr<const CLineModelResult> lmresult = std::dynamic_pointer_cast<const CLineModelResult>( linemodel_fpb.getResult() );
                    
        ChisquareArray chisquares = BuildChisquareArray(lmresult,
                                                        m_opt_rigidity,
                                                        m_opt_pdfcombination,
                                                        m_opt_stronglinesprior,
                                                        m_opt_haPrior,
                                                        m_opt_euclidNHaEmittersPriorStrength,
                                                        m_opt_modelZPriorStrength);
 

        //TODO deal with the case lmresult->Redshifts=1
        Int32 extremacount = 5;
        COperatorPdfz pdfz(m_opt_pdfcombination,
                            2*m_opt_secondpass_halfwindowsize, // peak separation
                            m_opt_candidatesLogprobaCutThreshold,
                            extremacount,
                            "FPB");

        std::shared_ptr<CPdfCandidateszResult> candResult = pdfz.Compute(chisquares, false);
        
        linemodel_fpb.SetFirstPassCandidates(candResult->m_ranked_candidates);

        //**************************************************
        //COMBINE CANDIDATES
        //**************************************************
        m_linemodel.Combine_firstpass_candidates(linemodel_fpb.m_firstpass_extremaResult);

    }

    //**************************************************
    //SECOND PASS
    //**************************************************
    bool skipSecondPass = (m_opt_skipsecondpass=="yes");
    if(!skipSecondPass)
    {
        Int32 retSecondPass = m_linemodel.ComputeSecondPass(spc,
                                                          tplCatalog,
                                                          tplCategoryList,
                                                          m_calibrationPath,
                                                          restraycatalog,
                                                          m_opt_linetypefilter,
                                                          m_opt_lineforcefilter,
                                                          lambdaRange,
                                                          m_opt_fittingmethod,
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
                                                          m_opt_abs_velocity_fit_step,
                                                          m_opt_secondpass_continuumfit);
        if( retSecondPass!=0 )
        {
            Log.LogError( "Linemodel, second pass failed. Aborting" );
            return false;
        }
    }else{
        m_linemodel.m_secondpass_parameters_extremaResult = *m_linemodel.m_firstpass_extremaResult;
    }


    //read it as constant to save it
    std::shared_ptr<const CLineModelResult> result = std::dynamic_pointer_cast<const CLineModelResult>( m_linemodel.getResult() );


    if( !result )
    {
        Log.LogError("%s: Failed to get linemodel result",__func__);
        throw runtime_error("Failed to get linemodel result");
    }else{
        //save linemodel chisquare results
        dataStore.StoreScopedGlobalResult( scopeStr.c_str(), result );

        //don't save linemodel extrema results, since will change with pdf computation

        //save linemodel firstpass extrema results
        std::string firstpassExtremaResultsStr=scopeStr.c_str();
        firstpassExtremaResultsStr.append("_firstpass_extrema");
        //Log.LogError("Linemodel, saving firstpass extrema results: %s", firstpassExtremaResultsStr.c_str());
        dataStore.StoreScopedGlobalResult( firstpassExtremaResultsStr.c_str(), m_linemodel.m_firstpass_extremaResult );

        //save linemodel firstpass extrema B results
        if(enableFirstpass_B)
        {
            std::string firstpassbExtremaResultsStr=scopeStr.c_str();
            firstpassbExtremaResultsStr.append("_firstpassb_extrema");
            //Log.LogError("Linemodel, saving firstpassb extrema results: %s", firstpassExtremaResultsStr.c_str());
            dataStore.StoreScopedGlobalResult( firstpassbExtremaResultsStr.c_str(), linemodel_fpb.m_firstpass_extremaResult );
        }

    }

    return true;
}
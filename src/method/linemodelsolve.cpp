#include <epic/redshift/method/linemodelsolve.h>

#include <epic/core/log/log.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/linemodel.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/datastore.h>

using namespace NSEpic;
using namespace std;

/** 
 * \brief Empty constructor.
 **/
CLineModelSolve::CLineModelSolve()
{

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
    desc.append("\tparam: linemodel.continuumcomponent = {""fromspectrum"", ""nocontinuum"", ""zero""}\n");
    desc.append("\tparam: linemodel.linewidthtype = {""psfinstrumentdriven"", ""zdriven"", ""fixedvelocity"", ""fixed""}\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.rules = {""all"", ""balmer"", ""strongweak"", ""oiiratio"", ""no""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");

    return desc;

}

/**
 * \brief Calls the Solve method and returns a new "result" object.
 * Call Solve.
 * Return a pointer to an empty CLineModelSolveResult. (The results for Linemodel will reside in the linemodel.linemodel result).
 **/
std::shared_ptr<const CLineModelSolveResult> CLineModelSolve::Compute( CDataStore& dataStore,
								       const CSpectrum& spc,
								       const CSpectrum& spcWithoutCont,
								       const CRayCatalog& restraycatalog,
								       const TFloat64Range& lambdaRange,
								       const TFloat64List& redshifts )
{
    CDataStore::CAutoScope resultScope( dataStore, "linemodelsolve" );
    Solve( dataStore, spc, spcWithoutCont, restraycatalog, lambdaRange, redshifts);
    return std::shared_ptr<const CLineModelSolveResult>( new CLineModelSolveResult() );
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
			     const CRayCatalog& restraycatalog,
                             const TFloat64Range& lambdaRange,
			     const TFloat64List& redshifts )
{
    std::string scopeStr = "linemodel";

    CSpectrum _spc = spc;
    CSpectrum _spcContinuum = spc;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract( spcWithoutCont.GetFluxAxis() );
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;


    std::string opt_linetypefilter;
    dataStore.GetScopedParam( "linemodel.linetypefilter", opt_linetypefilter, "no" );
    std::string opt_lineforcefilter;
    dataStore.GetScopedParam( "linemodel.lineforcefilter", opt_lineforcefilter, "S" );
    std::string opt_fittingmethod;
    dataStore.GetScopedParam( "linemodel.fittingmethod", opt_fittingmethod, "hybrid" );
    std::string opt_continuumcomponent;
    dataStore.GetScopedParam( "linemodel.continuumcomponent", opt_continuumcomponent, "fromspectrum" );
    std::string opt_lineWidthType;
    dataStore.GetScopedParam( "linemodel.linewidthtype", opt_lineWidthType, "fixedvelocity" );
    Float64 opt_resolution;
    dataStore.GetScopedParam( "linemodel.instrumentresolution", opt_resolution, 2350.0 );
    Float64 opt_velocity_emission;
    dataStore.GetScopedParam( "linemodel.velocityemission", opt_velocity_emission, 100.0 );
    Float64 opt_velocity_absorption;
    dataStore.GetScopedParam( "linemodel.velocityabsorption", opt_velocity_absorption, 300.0 );
    std::string opt_continuumreest;
    dataStore.GetScopedParam( "linemodel.continuumreestimation", opt_continuumreest, "no" );
    std::string opt_rules;
    dataStore.GetScopedParam( "linemodel.rules", opt_rules, "all" );
    Float64 opt_extremacount;
    dataStore.GetScopedParam( "linemodel.extremacount", opt_extremacount, 10.0 );

    Log.LogInfo( "Linemodel parameters:");
    Log.LogInfo( "    -linetypefilter: %s", opt_linetypefilter.c_str());
    Log.LogInfo( "    -lineforcefilter: %s", opt_lineforcefilter.c_str());
    Log.LogInfo( "    -fittingmethod: %s", opt_fittingmethod.c_str());
    Log.LogInfo( "    -linewidthtype: %s", opt_lineWidthType.c_str());
    Log.LogInfo( "    -rules: %s", opt_rules.c_str());
    if(opt_lineWidthType=="fixedvelocity"){
        Log.LogInfo( "    -instrumentresolution: %.2f", opt_resolution);
        Log.LogInfo( "    -velocity emission: %.2f", opt_velocity_emission);
        Log.LogInfo( "    -velocity absorption: %.2f", opt_velocity_absorption);
    }
    Log.LogInfo( "    -continuumreestimation: %s", opt_continuumreest.c_str());
    Log.LogInfo( "    -extremacount: %.3f", opt_extremacount);

    // Compute merit function
    COperatorLineModel linemodel;
    auto  result = linemodel.Compute( dataStore,
				      _spc,
				      _spcContinuum,
				      restraycatalog,
				      opt_linetypefilter,
				      opt_lineforcefilter,
				      lambdaRange,
				      redshifts,
				      opt_extremacount,
				      opt_fittingmethod,
				      opt_continuumcomponent,
				      opt_lineWidthType,
				      opt_resolution,
				      opt_velocity_emission,
				      opt_velocity_absorption,
				      opt_continuumreest,
				      opt_rules );

    if( !result )
    {
        //Log.LogInfo( "Failed to compute linemodel");
        return false;
    }else{
        // Store results
        dataStore.StoreScopedGlobalResult( scopeStr.c_str(), result );
    }

    return true;
}

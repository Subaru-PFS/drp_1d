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
    desc.append("\tparam: linemodel.linewidthtype = {""instrumentdriven"", ""combined"", ""nispsim2016"", ""fixed""}\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.rules = {""all"", ""balmer"", ""strongweak"", ""oiiratio"", ""ciiiratio"", ""no""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");
    desc.append("\tparam: linemodel.velocityfit = {""yes"", ""no""}\n");
    desc.append("\tparam: linemodel.fastfitlargegridstep = <float value>, deactivated if negative or zero\n");



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
    dataStore.GetScopedParam( "linemodel.fastfitlargegridstep", m_opt_twosteplargegridstep, 0.001 );
    dataStore.GetScopedParam( "linemodel.continuumcomponent", m_opt_continuumcomponent, "fromspectrum" );
    dataStore.GetScopedParam( "linemodel.linewidthtype", m_opt_lineWidthType, "combined" );
    dataStore.GetScopedParam( "linemodel.instrumentresolution", m_opt_resolution, 2350.0 );
    dataStore.GetScopedParam( "linemodel.velocityemission", m_opt_velocity_emission, 100.0 );
    dataStore.GetScopedParam( "linemodel.velocityabsorption", m_opt_velocity_absorption, 300.0 );
    dataStore.GetScopedParam( "linemodel.velocityfit", m_opt_velocityfit, "yes" );
    dataStore.GetScopedParam( "linemodel.continuumreestimation", m_opt_continuumreest, "no" );
    dataStore.GetScopedParam( "linemodel.rules", m_opt_rules, "all" );
    dataStore.GetScopedParam( "linemodel.extremacount", m_opt_extremacount, 10.0 );

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
    }else if(m_opt_lineWidthType=="nispsim2016"){
        Log.LogInfo( "    -velocity emission: %.2f", m_opt_velocity_emission);
        Log.LogInfo( "    -velocity absorption: %.2f", m_opt_velocity_absorption);
        Log.LogInfo( "    -velocity fit: %s", m_opt_velocityfit.c_str());
    }
    Log.LogInfo( "    -rules: %s", m_opt_rules.c_str());
    Log.LogInfo( "    -continuumreestimation: %s", m_opt_continuumreest.c_str());
    Log.LogInfo( "    -extremacount: %.3f", m_opt_extremacount);
    Log.LogInfo( "    -fastfitlargegridstep: %.6f", m_opt_twosteplargegridstep);

    return true;
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
    Bool storeResult = false;
    CDataStore::CAutoScope resultScope( dataStore, "linemodelsolve" );

    PopulateParameters( dataStore );
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



    // Compute merit function
    COperatorLineModel linemodel;
    auto  result = linemodel.Compute( dataStore,
				      _spc,
                      _spcContinuum,
				      restraycatalog,
                      m_opt_linetypefilter,
                      m_opt_lineforcefilter,
				      lambdaRange,
				      redshifts,
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
                      m_opt_twosteplargegridstep);

    if( !result )
    {
        //Log.LogInfo( "Failed to compute linemodel");
        return false;
    }else{
        // Store linemodel chisquare results
        dataStore.StoreScopedGlobalResult( scopeStr.c_str(), result );
        //save linemodel fitting and spectrum-model results
        linemodel.storeGlobalModelResults(dataStore);
    }

    return true;
}

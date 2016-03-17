#include <epic/redshift/method/linemodeltplshapesolve.h>

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
CLineModelTplshapeSolve::CLineModelTplshapeSolve()
{

}

/**
 * \brief Empty destructor.
 **/
CLineModelTplshapeSolve::~CLineModelTplshapeSolve()
{

}

/**
 * \brief Returns a string describing the names and allowed values for the parameters of the Linemodel method.
 **/
const std::string CLineModelTplshapeSolve::GetDescription()
{
    std::string desc;

    desc = "Method linemodeltplshape:\n";

    desc.append("\tparam: linemodel.linetypefilter = {""no"", ""E"", ""A""}\n");
    desc.append("\tparam: linemodel.lineforcefilter = {""no"", ""S""}\n");
    desc.append("\tparam: linemodel.fittingmethod = {""hybrid"", ""individual""}\n");
    desc.append("\tparam: linemodel.continuumcomponent = {""fromspectrum"", ""nocontinuum"", ""zero""}\n");
    desc.append("\tparam: linemodel.linewidthtype = {""instrumentdriven"", ""combined"", ""fixed""}\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.rules = {""all"", ""balmer"", ""strongweak"", ""oiiratio"", ""ciiiratio"", ""no""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");
    desc.append("\tparam: linemodel.velocityfit = {""yes"", ""no""}\n");



    return desc;

}

/**
 * \brief
 * Populates the method parameters from the dataStore into the class members
 * Returns true if successful, false otherwise
 **/
Bool CLineModelTplshapeSolve::PopulateParameters( CDataStore& dataStore )
{
    dataStore.GetScopedParam( "linemodel.linetypefilter", m_opt_linetypefilter, "no" );
    dataStore.GetScopedParam( "linemodel.lineforcefilter", m_opt_lineforcefilter, "S" );
    dataStore.GetScopedParam( "linemodel.fittingmethod", m_opt_fittingmethod, "hybrid" );
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
    }
    Log.LogInfo( "    -rules: %s", m_opt_rules.c_str());
    Log.LogInfo( "    -continuumreestimation: %s", m_opt_continuumreest.c_str());
    Log.LogInfo( "    -extremacount: %.3f", m_opt_extremacount);

    return true;
}

/**
 * \brief Calls the Solve method and returns a new "result" object.
 * Call Solve.
 * Return a pointer to an empty CLineModelTplshapeSolveResult. (The results for Linemodel will reside in the linemodel.linemodel result).
 **/
std::shared_ptr<const CLineModelTplshapeSolveResult> CLineModelTplshapeSolve::Compute( CDataStore& dataStore,
								       const CSpectrum& spc,
								       const CSpectrum& spcWithoutCont,
								       const CRayCatalog& restraycatalog,
								       const TFloat64Range& lambdaRange,
								       const TFloat64List& redshifts )
{
    Bool storeResult = false;
    CDataStore::CAutoScope resultScope( dataStore, "linemodeltplshapesolve" );

    PopulateParameters( dataStore );
    Solve( dataStore, spc, spcWithoutCont, restraycatalog, lambdaRange, redshifts);
    return std::shared_ptr<const CLineModelTplshapeSolveResult>( new CLineModelTplshapeSolveResult() );
}

/**
 * \brief
 * Create a continuum object by subtracting spcWithoutCont from the spc.
 * Configure the opt_XXX variables from the dataStore scope parameters.
 * LogInfo the opt_XXX values.
 * Create a COperatorLineModel, call its Compute method. 
 * If that returned true, store results.
 **/
Bool CLineModelTplshapeSolve::Solve( CDataStore& dataStore,
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

    //load the catalogs list from the files in the tplshape-catalogs folder : tplshapeCatalogDir
    namespace fs = boost::filesystem;
    fs::path tplshapeCatalogDir("/home/aschmitt/data/vuds/VUDS_flag3_4/amazed/linecatalogs/linecatalogs_tplshape_ExtendedGalaxyEL2_20160317");
    fs::directory_iterator end_iter;
    std::vector<std::string> tplshapeCatalogList;
    if ( fs::exists(tplshapeCatalogDir) && fs::is_directory(tplshapeCatalogDir))
    {
      for( fs::directory_iterator dir_iter(tplshapeCatalogDir) ; dir_iter != end_iter ; ++dir_iter)
      {
        if (fs::is_regular_file(dir_iter->status()) )
        {
          tplshapeCatalogList.push_back(dir_iter->path().c_str());
        }
      }
    }
    Log.LogInfo( "Linemodel tplshape - Found %d tplshaped catalogs", tplshapeCatalogList.size());

    //prepare the combined merit function
    std::shared_ptr<CLineModelResult> resultCombined = std::shared_ptr<CLineModelResult>( new CLineModelResult() );

    // Compute merit function for each tplshape_catalog
    for(Int32 k=0; k<tplshapeCatalogList.size(); k++)
    {
        //get line catalog
        CRayCatalog lineCatalog;
        Bool rValue = lineCatalog.Load( tplshapeCatalogList[k].c_str() );
        if( !rValue )
        {
            Log.LogInfo( "Failed to load tplshape catalog: %s", tplshapeCatalogList[k].c_str());
            continue;
        }
        else
        {
            Log.LogInfo( "Loaded tplshape: %s", tplshapeCatalogList[k].c_str());
        }
        COperatorLineModel linemodel;
        auto  result = linemodel.Compute( dataStore,
                                          _spc,
                                          _spcContinuum,
                                          lineCatalog,
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
                                          m_opt_velocityfit);

        if( !result )
        {
            //Log.LogInfo( "Failed to compute linemodel");
            return false;
        }else{
            auto resultPtr = std::dynamic_pointer_cast<const CLineModelResult>( result );

            if(resultCombined->Redshifts.size()<1)
            {
                Int32 nz = resultPtr->Redshifts.size();
                resultCombined->ChiSquare.resize( nz );
                resultCombined->Redshifts.resize( nz );
                resultCombined->LineModelSolutions.resize( nz );

                for( Int32 i=0; i<nz; i++ )
                {
                    resultCombined->ChiSquare[i] = resultPtr->ChiSquare[i];
                    resultCombined->Redshifts[i] = resultPtr->Redshifts[i];
                    resultCombined->LineModelSolutions[i] = resultPtr->LineModelSolutions[i];
                }

                resultCombined->ResizeExtremaResults(1);
                resultCombined->Extrema[0] = resultPtr->Extrema[0];

            }
            else
            {

                Int32 n = resultCombined->Redshifts.size();
                for( Int32 i=0; i<n; i++ )
                {
                    if( resultPtr->ChiSquare[i] < resultCombined->ChiSquare[i] )
                    {
                        resultCombined->ChiSquare[i] = resultPtr->ChiSquare[i];
                        resultCombined->LineModelSolutions[i] = resultPtr->LineModelSolutions[i];
                    }
                }

                //deal with the first extrema only, for now...
                Int32 extremaIdx = 0;
                Int32 solutionCombinedIdx=0;
                for ( UInt32 i2=0; i2<resultCombined->LineModelSolutions.size(); i2++)
                {
                    if( resultCombined->Redshifts[i2]==resultCombined->Extrema[extremaIdx] )
                    {
                        solutionCombinedIdx = i2;
                        break;
                    }
                }
                Int32 solutionIdx=0;
                for ( UInt32 i2=0; i2<resultPtr->LineModelSolutions.size(); i2++)
                {
                    if( resultPtr->Redshifts[i2]==resultPtr->Extrema[extremaIdx] )
                    {
                        solutionIdx = i2;
                        break;
                    }
                }

                if( resultPtr->ChiSquare[solutionIdx] < resultCombined->ChiSquare[solutionCombinedIdx] )
                {
                    resultCombined->Extrema[extremaIdx] = resultPtr->Extrema[extremaIdx];
                }
            }
        }
    }

    // Store results
    dataStore.StoreScopedGlobalResult( scopeStr.c_str(), resultCombined );

    return true;
}

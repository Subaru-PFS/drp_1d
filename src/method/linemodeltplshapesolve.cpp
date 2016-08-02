#include <epic/redshift/method/linemodeltplshapesolve.h>

#include <epic/core/log/log.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/linemodel.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/datastore.h>

#include <fstream>

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
                                       const CTemplateCatalog& tplCatalog,
                                       const TStringList& tplCategoryList,
								       const CRayCatalog& restraycatalog,
								       const TFloat64Range& lambdaRange,
								       const TFloat64List& redshifts )
{
    Bool storeResult = false;
    CDataStore::CAutoScope resultScope( dataStore, "linemodeltplshapesolve" );

    PopulateParameters( dataStore );


    CSpectrum _spc = spc;
    CSpectrum _spcContinuum = spc;
    CSpectrumFluxAxis spcfluxAxis = _spcContinuum.GetFluxAxis();
    spcfluxAxis.Subtract( spcWithoutCont.GetFluxAxis() );
    CSpectrumFluxAxis& sfluxAxisPtr = _spcContinuum.GetFluxAxis();
    sfluxAxisPtr = spcfluxAxis;


    //load the catalogs list from the files in the tplshape-catalogs folder : tplshapeCatalogDir
    namespace fs = boost::filesystem;
    std::string dirPath = "/home/aschmitt/data/vuds/VUDS_flag3_4/amazed/linecatalogs/linecatalogs_tplshape_ExtendedTemplatesMarch2016_v2_20160731_B10H";
    fs::path tplshapeCatalogDir(dirPath.c_str());

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

    //load the velocities list for all the catalogs
    fs::path tplshapeVelocitiesDir = tplshapeCatalogDir/"velocities/";
    std::vector<std::string> tplshapeVelocitiesList;
    if ( fs::exists(tplshapeVelocitiesDir) && fs::is_directory(tplshapeVelocitiesDir))
    {
      for( fs::directory_iterator dir_iter(tplshapeVelocitiesDir) ; dir_iter != end_iter ; ++dir_iter)
      {
        if (fs::is_regular_file(dir_iter->status()) )
        {
          tplshapeVelocitiesList.push_back(dir_iter->path().c_str());
        }
      }
    }
    Log.LogInfo( "Linemodel tplshape - Found %d tplshaped velocities files", tplshapeVelocitiesList.size());

    //loop on the templates
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );

            //find the linecatalog-tplshaped corresponding to the template name
            Int32 kctlg = -1;
            for(Int32 k=0; k<tplshapeCatalogList.size(); k++)
            {
                std::string tplname = tpl.GetName();
                std::string ctlgname = tplshapeCatalogList[k];
                std::size_t foundstra = ctlgname.find(tplname.c_str());
                if (foundstra==std::string::npos){
                    continue;
                }
                kctlg = k;
            }

            if(kctlg<0)
            {
                Log.LogError( "Failed to match tpl with tplshape catalog: %s", tpl.GetName().c_str());
                continue;
            }
            //get line catalog
            CRayCatalog lineCatalog;
            Bool rValue = lineCatalog.Load( tplshapeCatalogList[kctlg].c_str() );
            if( !rValue )
            {
                Log.LogError( "Failed to load tplshape catalog: %s", tplshapeCatalogList[kctlg].c_str());
                continue;
            }
            else
            {
                Log.LogInfo( "Loaded tplshape: %s", tplshapeCatalogList[kctlg].c_str());
            }

            if(0)
            {
                //find the velocities-tplshaped corresponding to the template name
                Int32 kvel = -1;
                for(Int32 k=0; k<tplshapeVelocitiesList.size(); k++)
                {
                    std::string tplname = tpl.GetName();
                    std::string velname = tplshapeVelocitiesList[k];
                    std::size_t foundstra = velname.find(tplname.c_str());
                    if (foundstra==std::string::npos){
                        continue;
                    }
                    kvel = k;
                }

                if(kvel<0)
                {
                    Log.LogError( "Failed to match tpl with tplshape velocities: %s", tpl.GetName().c_str());
                    continue;
                }
                //get velocities from file
                Float64 elv=100.0;
                Float64 alv=300.0;
                bool ret = LoadVelocities(tplshapeVelocitiesList[kvel].c_str(), elv, alv);
                if( !ret )
                {
                    Log.LogError( "Failed to load tplshape velocities: %s", tplshapeVelocitiesList[kvel].c_str());
                    continue;
                }
                else
                {
                    Log.LogInfo( "Loaded tplshape velocities: %s", tplshapeVelocitiesList[kvel].c_str());
                    m_opt_velocity_emission = elv;
                    m_opt_velocity_absorption = alv;
                }
            }

            Solve( dataStore, _spc, _spcContinuum, tpl, lineCatalog, lambdaRange, redshifts);
            storeResult = true;
        }
    }

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
                 const CSpectrum& spcCont,
                 const CTemplate& tpl,
                 const CRayCatalog& lineCatalog,
                 const TFloat64Range& lambdaRange,
			     const TFloat64List& redshifts )
{
    std::string scopeStr = "linemodel";

    COperatorLineModel linemodel;
    auto  result = linemodel.Compute( dataStore,
                                      spc,
                                      spcCont,
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
        Log.LogInfo( "Failed to compute linemodeltplshape");
        return false;
    }else{
        //save linemodel chisquare results
        dataStore.StoreScopedPerTemplateResult( tpl, scopeStr.c_str(), result );
        //save linemodel fitting and spectrum-model results
        linemodel.storePerTemplateModelResults(dataStore, tpl);
    }


    return true;
}

Bool CLineModelTplshapeSolve::LoadVelocities( const char* filePath, Float64& elv, Float64& alv )
{
    ifstream file;

    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
        return false;

    string line;

    // Read file line by line
    Int32 readNums = 0;
    while( getline( file, line ) )
    {
        if(readNums==0)
        {
            elv = std::stod(line);
        }else if(readNums==1)
        {
            alv = std::stod(line);
        }
        readNums++;
    }

    if(readNums!=2)
    {
        return false;
    }

    return true;
}

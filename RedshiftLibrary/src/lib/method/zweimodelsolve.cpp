#include <RedshiftLibrary/method/zweimodelsolve.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/noise/flat.h>
#include <RedshiftLibrary/noise/fromfile.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

#include <RedshiftLibrary/statistics/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/operator/pdfLogresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iostream>

namespace bfs = boost::filesystem;

using namespace NSEpic;
using namespace std;
using namespace boost;


/** 
 * \brief Empty constructor.
 **/
CZweiModelSolve::CZweiModelSolve(string calibrationPath)
{
    m_calibrationPath = calibrationPath;
}

/**
 * \brief Empty destructor.
 **/
CZweiModelSolve::~CZweiModelSolve()
{

}

/**
 * \brief Returns a string describing the names and allowed values for the parameters of the Linemodel method.
 **/
const std::string CZweiModelSolve::GetDescription()
{
    std::string desc;

    desc = "Method zweimodel:\n";

    desc.append("\tparam: zweimodel.linetypefilter = {""no"", ""E"", ""A""}\n");
    desc.append("\tparam: zweimodel.lineforcefilter = {""no"", ""S""}\n");
    desc.append("\tparam: zweimodel.fittingmethod = {""hybrid"", ""individual""}\n");
    desc.append("\tparam: zweimodel.continuumcomponent = {""fromspectrum"", ""tplfit"", ""nocontinuum"", ""zero""}\n");
    desc.append("\tparam: zweimodel.rigidity = {""rules"", ""tplcorr"", ""tplshape""}\n");
    desc.append("\tparam: zweimodel.linewidthtype = {""instrumentdriven"", ""velocitydriven"",  ""combined"",  ""nispvsspsf201707"", ""fixed""}\n");
    desc.append("\tparam: zweimodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: zweimodel.velocityemission = <float value>\n");
    desc.append("\tparam: zweimodel.velocityabsorption = <float value>\n");
    desc.append("\tparam: zweimodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: zweimodel.rules = {""all"", ""balmer"", ""strongweak"", ""superstrong"", ""ratiorange"", ""ciiiratio"", ""no""}\n");
    desc.append("\tparam: zweimodel.extremacount = <float value>\n");
    desc.append("\tparam: zweimodel.velocityfit = {""yes"", ""no""}\n");
    desc.append("\tparam: zweimodel.emvelocityfitmin = <float value>\n");
    desc.append("\tparam: zweimodel.emvelocityfitmax = <float value>\n");
    desc.append("\tparam: zweimodel.absvelocityfitmin = <float value>\n");
    desc.append("\tparam: zweimodel.absvelocityfitmax = <float value>\n");
    desc.append("\tparam: zweimodel.firstpass.largegridstep = <float value>, deactivated if negative or zero\n");
    desc.append("\tparam: zweimodel.pdfcombination = {""marg"", ""bestchi2""}\n");
    desc.append("\tparam: zweimodel.stronglinesprior = <float value>, penalization factor = positive value or -1 to deactivate\n");
    desc.append("\tparam: zweimodel.saveintermediateresults = {""yes"", ""no""}\n");



    return desc;

}

/**
 * \brief
 * Populates the method parameters from the dataStore into the class members
 * Returns true if successful, false otherwise
 **/
Bool CZweiModelSolve::PopulateParameters( CDataStore& dataStore )
{
    dataStore.GetScopedParam( "zweimodel.linetypefilter", m_opt_linetypefilter, "no" );
    dataStore.GetScopedParam( "zweimodel.lineforcefilter", m_opt_lineforcefilter, "no" );
    dataStore.GetScopedParam( "zweimodel.fittingmethod", m_opt_fittingmethod, "hybrid" );
    dataStore.GetScopedParam( "zweimodel.firstpass.largegridstep", m_opt_firstpass_largegridstep, 0.001 );
    std::string redshiftSampling;
    dataStore.GetParam( "redshiftsampling", redshiftSampling, "lin" ); //TODO: sampling in log cannot be used for now as zqual descriptors assume constant dz.
    if(redshiftSampling=="log")
    {
        m_opt_firstpass_largegridsampling = "log";
    }else{
        m_opt_firstpass_largegridsampling = "lin";
    }
    Log.LogDetail( "    firstpass - largegridsampling (auto set from redshiftsampling param.): %s", m_opt_firstpass_largegridsampling.c_str());

    dataStore.GetScopedParam( "zweimodel.continuumcomponent", m_opt_continuumcomponent, "fromspectrum" );
    dataStore.GetScopedParam( "zweimodel.rigidity", m_opt_rigidity, "rules" );
    if(m_opt_rigidity=="tplshape")
    {
        dataStore.GetScopedParam( "zweimodel.tplratio_catalog", m_opt_tplratio_reldirpath, "linecatalogs_tplshapes/linecatalogs_tplshape_ExtendedTemplatesJan2017v3_20170602_B14C_v3_emission" );
        dataStore.GetScopedParam( "zweimodel.tplratio_ismfit", m_opt_tplratio_ismfit, "yes" );
    }
    dataStore.GetScopedParam( "zweimodel.offsets_catalog", m_opt_offsets_reldirpath, "linecatalogs_offsets/offsetsCatalogs_20170410_m150" );

    dataStore.GetScopedParam( "zweimodel.linewidthtype", m_opt_lineWidthType, "velocitydriven" );
    dataStore.GetScopedParam( "zweimodel.instrumentresolution", m_opt_resolution, 2350.0 );
    dataStore.GetScopedParam( "zweimodel.velocityemission", m_opt_velocity_emission, 100.0 );
    dataStore.GetScopedParam( "zweimodel.velocityabsorption", m_opt_velocity_absorption, 300.0 );
    dataStore.GetScopedParam( "zweimodel.velocityfit", m_opt_velocityfit, "yes" );
    if(m_opt_velocityfit=="yes"){
        dataStore.GetScopedParam( "zweimodel.emvelocityfitmin", m_opt_em_velocity_fit_min, 20.0 );
        dataStore.GetScopedParam( "zweimodel.emvelocityfitmax", m_opt_em_velocity_fit_max, 500.0 );
        dataStore.GetScopedParam( "zweimodel.absvelocityfitmin", m_opt_abs_velocity_fit_min, 150.0 );
        dataStore.GetScopedParam( "zweimodel.absvelocityfitmax", m_opt_abs_velocity_fit_max, 500.0 );
    }
    dataStore.GetScopedParam( "zweimodel.continuumreestimation", m_opt_continuumreest, "no" );
    dataStore.GetScopedParam( "zweimodel.rules", m_opt_rules, "all" );
    dataStore.GetScopedParam( "zweimodel.extremacount", m_opt_extremacount, 10.0 );
    dataStore.GetScopedParam( "zweimodel.stronglinesprior", m_opt_stronglinesprior, 1e-1);
    dataStore.GetScopedParam( "zweimodel.pdfcombination", m_opt_pdfcombination, "marg");
    dataStore.GetScopedParam( "zweimodel.saveintermediateresults", m_opt_saveintermediateresults, "no");

    //Auto-correct fitting method
    //std::string forcefittingmethod = "ones";
    std::string forcefittingmethod = "individual";
    if(m_opt_rigidity=="tplshape" && m_opt_fittingmethod != forcefittingmethod)
    {
        m_opt_fittingmethod = forcefittingmethod;
        dataStore.SetScopedParam("zweimodel.fittingmethod", m_opt_fittingmethod);
        Log.LogInfo( "LineModel fitting method auto-correct due to tplshape rigidity");

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
        Log.LogInfo( "    -abs velocity fit min : %.1f", m_opt_abs_velocity_fit_min);
        Log.LogInfo( "    -abs velocity fit max : %.1f", m_opt_abs_velocity_fit_max);
    }

    Log.LogInfo( "    -rigidity: %s", m_opt_rigidity.c_str());
    if(m_opt_rigidity=="rules"){
        Log.LogInfo( "    -rules: %s", m_opt_rules.c_str());
    }else if(m_opt_rigidity=="tplshape")
    {
        Log.LogInfo( "      -tplratio_catalog", m_opt_tplratio_reldirpath.c_str());
        Log.LogInfo( "      -tplratio_ismfit", m_opt_tplratio_ismfit.c_str());
    }
    Log.LogInfo( "      -offsets_catalog", m_opt_offsets_reldirpath.c_str());

    Log.LogInfo( "    -continuumcomponent: %s", m_opt_continuumcomponent.c_str());
    Log.LogInfo( "    -continuumreestimation: %s", m_opt_continuumreest.c_str());
    Log.LogInfo( "    -extremacount: %.3f", m_opt_extremacount);
    Log.LogInfo( "    -first pass:");
    Log.LogInfo( "      -largegridstep: %.6f", m_opt_firstpass_largegridstep);

    Log.LogInfo( "    -pdf-stronglinesprior: %e", m_opt_stronglinesprior);
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
 * Return a pointer to an empty CLineModelSolveResult. (The results for Linemodel will reside in the linemodelsolve.zweimodel result).
 * NB: keep the linemodelsolve scope for now, so that the linemodelsolveresult can be used
 **/
std::shared_ptr<CLineModelSolveResult> CZweiModelSolve::Compute( CDataStore& dataStore,
								       const CSpectrum& spc,
								       const CSpectrum& spcWithoutCont,
                                       const CTemplateCatalog& tplCatalog,
                                       const TStringList& tplCategoryList,
								       const CRayCatalog& restraycatalog,
								       const TFloat64Range& lambdaRange,
								       const TFloat64List& redshifts )
{
    CDataStore::CAutoScope resultScope( dataStore, "linemodelsolve" );

    PopulateParameters( dataStore );
    Int32 retSolve = Solve( dataStore, spc, spcWithoutCont, tplCatalog, tplCategoryList, restraycatalog, lambdaRange, redshifts);

    if(retSolve){

        /* ------------------------  COMPUTE POSTMARG PDF  --------------------------  */
        Log.LogInfo("zweimodelsolve: Pdfz computation");

        std::string scope = "linemodelsolve.linemodel";
        auto results = dataStore.GetGlobalResult( scope.c_str() );
        if(results.expired())
        {
            Log.LogError("zweimodelsolve: Unable to retrieve linemodel results");
            return NULL;
        }
        std::shared_ptr<const CLineModelResult> result = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );

        CombinePDF(dataStore, result, m_opt_rigidity, m_opt_pdfcombination, m_opt_stronglinesprior);

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
    }

    return std::shared_ptr<CLineModelSolveResult>( new CLineModelSolveResult() );
}


Int32 CZweiModelSolve::CombinePDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result, std::string opt_rigidity, std::string opt_combine, Float64 opt_stronglinesprior)
{
    bool zPriorStrongLinePresence = (opt_stronglinesprior>0.0);
    if(zPriorStrongLinePresence)
    {
        Log.LogInfo("Linemodel: Pdfz computation: StrongLinePresence prior enabled: factor=%e", opt_stronglinesprior);
    }else{
        Log.LogInfo("Linemodel: Pdfz computation: StrongLinePresence prior disabled");
    }

    Log.LogInfo("Linemodel: Pdfz computation");
    std::shared_ptr<CPdfMargZLogResult> postmargZResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
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
        std::shared_ptr<CPdfLogResult> zPrior = std::shared_ptr<CPdfLogResult>(new CPdfLogResult());
        zPrior->SetSize(result->Redshifts.size());
        for ( UInt32 k=0; k<result->Redshifts.size(); k++)
        {
            zPrior->Redshifts[k] = result->Redshifts[k];
        }
        if(zPriorStrongLinePresence)
        {
            UInt32 lineTypeFilter = 1;// for emission lines only
            std::vector<bool> strongLinePresence = result->GetStrongLinesPresence(lineTypeFilter, result->LineModelSolutions);
            zPrior->valProbaLog = pdfz.GetStrongLinePresenceLogZPrior(strongLinePresence, opt_stronglinesprior);
        }else{
            zPrior->valProbaLog = pdfz.GetConstantLogZPrior(result->Redshifts.size());
        }

        //correct chi2 if necessary: todo add switch
        TFloat64List logLikelihoodCorrected(result->ChiSquare.size(), DBL_MAX);
        for ( UInt32 k=0; k<result->Redshifts.size(); k++)
        {
            logLikelihoodCorrected[k] = result->ChiSquare[k];// + result->ScaleMargCorrection[k];
        }
        retPdfz = pdfz.Compute(logLikelihoodCorrected, result->Redshifts, cstLog, zPrior->valProbaLog, logProba, logEvidence);
        if(retPdfz==0){
            store.StoreGlobalResult( "zPDF/logprior.logP_Z_data", zPrior);

            postmargZResult->countTPL = result->Redshifts.size(); // assumed 1 model per z
            postmargZResult->Redshifts.resize(result->Redshifts.size());
            postmargZResult->valProbaLog.resize(result->Redshifts.size());
            for ( UInt32 k=0; k<result->Redshifts.size(); k++)
            {
                postmargZResult->Redshifts[k] = result->Redshifts[k] ;
                postmargZResult->valProbaLog[k] = logProba[k];
            }
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


    if(retPdfz!=0)
    {
        Log.LogError("Linemodel: Pdfz computation failed");
    }else{
        store.StoreGlobalResult( "zPDF/logposterior.logMargP_Z_data", postmargZResult); //need to store this pdf with this exact same name so that zqual can load it. see zqual.cpp/ExtractFeaturesPDF
    }

    return 0;
}


Int32 CZweiModelSolve::SaveContinuumPDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result)
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
Int32 CZweiModelSolve::getVelocitiesFromRefFile( const char* filePath, std::string spcid, Float64& elv, Float64& alv )
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
 * Retrieve the true-redshift from a hardcoded ref file path
 * nb: this is a hack for development purposes/line measurement at zref
 * reverseInclusion=0 (default): spcId is searched to be included in the Ref-File-Id
 * reverseInclusion=1 : Ref-File-Id is searched to be included in the spcId
 **/
Int32 CZweiModelSolve::getValueFromRefFile( const char* filePath, std::string spcid, Int32 colID, Float64& zref, Int32 reverseInclusion )
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


/**
 * \brief
 * Create a continuum object by subtracting spcWithoutCont from the spc.
 * Configure the opt_XXX variables from the dataStore scope parameters.
 * LogInfo the opt_XXX values.
 * Create a COperatorLineModel, call its Compute method. 
 * If that returned true, store results.
 **/
Bool CZweiModelSolve::Solve( CDataStore& dataStore,
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

    // ---------------------------------------------------
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

    // ---------------------------------------------------
    //Hack: load the zref values
    std::vector<Float64> _redshifts;
    if(false)
    {
        Log.LogInfo( "Linemodel - hacking zref enabled");
        Float64 zref = -1.0;
        namespace fs = boost::filesystem;
        Int32 reverseInclusionForIdMatching = 0;
        bool computeOnZrange;
        /*
        // Euclid case for process at zref
        fs::path refFilePath("/home/aschmitt/amazed_cluster/datasets/euclid/euclidsim2016/EUC-TEST-TUGALSPC-2016-03_export20170302_3z4mag4lfhabins_TrueFlux/reference_correctedFastSim.list");
        //fs::path refFilePath("/sps/euclid/Users/schmitt/amazed_cluster/datasets/euclid/euclidsim2016/EUC-TEST-TUGALSPC-2016-03_export20170302_3z4mag4lfhabins_TrueFlux/reference_correctedFastSim.list");
        Int32 substring_start = 0;
        Int32 substring_n = spc.GetName().size();
        Int32 colId = 2; //starts at 1, so that id_column is usually 1
        reverseInclusionForIdMatching = 1;
        computeOnZrange = true;
        //*/
        //*
        // SDSS case for process at zref
        //fs::path refFilePath("/home/aschmitt/amazed_cluster/datasets/sdss/sdss_201707/reference_SDSS_spectra_bg10k.txt");
        //fs::path refFilePath("/sps/euclid/Users/schmitt/amazed_cluster/datasets/sdss/sdss_201707/reference_SDSS_spectra_bg10k.txt");
        fs::path refFilePath("/home/aschmitt/data/sdss/sdss_201707/reference_SDSS_spectra_bg10k.txt");
        Int32 substring_start = 5;
        Int32 substring_n = 15;
        Int32 colId = 2; //starts at 1, so that id_column is usually 1
        computeOnZrange = false;
        //*/
        /*
        // PFS case for process at zref
        fs::path refFilePath("/home/aschmitt/amazed_cluster/datasets/pfs/pfs6b2_201704_sim10k/pfs6b2_reference_20170523_filtSim3h.txt");
        Int32 substring_start = 0;
        Int32 substring_n = 18;
        std::string strTag = "wlines";
        std::size_t foundstra = spc.GetName().find(strTag.c_str());
        if (foundstra!=std::string::npos){
            substring_n = (Int32)foundstra;
        }else{
            Log.LogWarning( "Linemodel - hack - unable to find strTag=%s", strTag.c_str());
            return false;
        }

        Int32 colId = 2; //starts at 1, so that id_column is usually 1
        computeOnZrange = false;
        //*/

        if ( fs::exists(refFilePath) )
        {
            std::string spcSubStringId = spc.GetName().substr(substring_start, substring_n);
            Log.LogInfo( "Linemodel - hack - using substring %s", spcSubStringId.c_str());
            getValueFromRefFile( refFilePath.c_str(), spcSubStringId, colId, zref, reverseInclusionForIdMatching);
        }
        if(zref==-1)
        {
            Log.LogWarning( "Linemodel - hack - unable to find zref!");
            return false;
        }

        if(computeOnZrange) //computing only on zref, or on a zrange around zref
        {
            Float64 deltaZrangeHalf = 0.5e-2; //override zrange
            Float64 stepZ = 1e-5;
            Float64 nStepsZ = deltaZrangeHalf*2/stepZ+1;
            for(Int32 kz=0; kz<nStepsZ; kz++)
            {
                Float64 _z = zref + kz*stepZ - deltaZrangeHalf;
                _redshifts.push_back(_z);
            }
            Log.LogInfo( "Linemodel - hack - zmin=%.5f, zmax=%.5f, zstep=%.5f", _redshifts[0], _redshifts[_redshifts.size()-1], stepZ);
        }else{
            _redshifts.push_back(zref);
        }
        m_opt_extremacount = 1; //override nextrema count
        m_opt_firstpass_largegridstep = 0; //override fastfitlargegrid
        Log.LogInfo( "Linemodel - hack - Loaded zref for spc %s : zref=%f", spc.GetName().c_str(), zref);
    }else{
        _redshifts = redshifts;
    }


    // ---------------------------------------------------
    //get the contaminant spectra/noise paths
    Log.LogInfo("===============================================");
    std::string spcContFileName="";
    std::string errorContFileName="";
    Float64 offsetLambdaContaminant=0.0;
    std::string contMatchFilePath = "/home/aschmitt/data/euclid/simulation2017-SC3_test_zweiroll/simulation2017-SC3_test_zweiroll_source1source3/simulation2017-SC3_test_zweiroll_source1_rolls/simulation2017-SC3_test_zweiroll_source1_source3_match.txt";
    Int32 reverseInclusionForIdMatching = 0;
    namespace fs = boost::filesystem;
    if ( fs::exists(contMatchFilePath) )
    {
        Int32 substring_start = 0;
        Int32 substring_n = spc.GetName().size();
        std::string strTag = "_F";
        std::size_t foundstra = spc.GetName().find(strTag.c_str());
        if (foundstra!=std::string::npos){
            substring_n = (Int32)foundstra;
        }else{
            Log.LogWarning( "Zweimodel - contaminant file search - unable to find strTag=%s", strTag.c_str());
            return false;
        }


        std::string spcSubStringId = spc.GetName().substr(substring_start, substring_n);
        Log.LogInfo( "Zweimodel - contaminant file search - using substring %s", spcSubStringId.c_str());
        Int32 colId=1;
        Bool retContFiles = getContaminantNameFromFile( contMatchFilePath.c_str(), spcSubStringId, colId, spcContFileName, errorContFileName, offsetLambdaContaminant, reverseInclusionForIdMatching);
        if(retContFiles==false)
        {
            Log.LogError( "Zweimodel - contaminant - failed to read cont. match file=%s", contMatchFilePath.c_str());
        }
        Log.LogInfo( "Zweimodel - contaminant - using cont spc : %s", spcContFileName.c_str());
        Log.LogInfo( "Zweimodel - contaminant - using cont error : %s", errorContFileName.c_str());
        Log.LogInfo( "Zweimodel - contaminant - using offsetLambdaContaminant_s2tos1 : %.2f", offsetLambdaContaminant);

    }else{
        Log.LogError( "Zweimodel - contaminant - unable to find cont. match file=%s", contMatchFilePath.c_str());
        return false;
    }

    fs::path spcdirPath = bfs::path(spc.GetFullPath()).branch_path();
    Log.LogInfo( "Zweimodel - contaminant - using spcdir : %s", spcdirPath.string().c_str());
    std::string spcContFilePath;
    std::string errorContFilePath;
    spcContFilePath = ( spcdirPath/bfs::path(spcContFileName) ).string();
    errorContFilePath = ( spcdirPath/bfs::path(errorContFileName) ).string();

    //load the contaminant data
    std::shared_ptr<CSpectrum> contSpectrum = std::shared_ptr<CSpectrum>( new CSpectrum() );
    std::string spcName="";
    if(spcContFilePath.c_str()!=NULL)
    {
        spcName = bfs::path( spcContFilePath ).stem().string() ;
    }
    std::string noiseName="";
    if(errorContFilePath.c_str()!=NULL)
    {
        noiseName = bfs::path( errorContFilePath ).stem().string() ;
    }
    contSpectrum->SetName(spcName.c_str());
    contSpectrum->SetFullPath(bfs::path( spcContFilePath ).string().c_str() );

    CSpectrumIOGenericReader reader;

    reader.Read( spcContFilePath.c_str(), *contSpectrum );
    Log.LogInfo("Zweimodel - contaminant - Successfully loaded contaminant spectrum file: (%s)", spcName.c_str() );

    // add noise if any or add flat noise
    if( errorContFilePath.c_str() == NULL )
    {
        CNoiseFlat noise;
        noise.SetStatErrorLevel( 1.0 );
        noise.AddNoise( *contSpectrum );
    }
    else
    {
        CNoiseFromFile noise;
        noise.SetNoiseFilePath( errorContFilePath.c_str(), reader );
        {
            Log.LogError("Zweimodel - contaminant - Failed to load noise spectrum");
            return false;
        }

        noise.AddNoise( *contSpectrum );
	Log.LogInfo("Zweimodel - contaminant - Successfully loaded input noise file:    (%s)", noiseName.c_str() );
    }
    Log.LogInfo("===============================================");


    // ---------------------------------------------------
    // Compute the S1 first pass without contaminant
    Log.LogInfo("Zweimodel - Computing first linemodel : S1");
    Log.LogInfo("Zweimodel - ===============================================");
    COperatorLineModel linemodel_s1;

    linemodel_s1.m_secondPass_velfit_dzInfLim = 4e-4;
    linemodel_s1.m_secondPass_velfit_dzSupLim = 4e-4;
    linemodel_s1.m_secondPass_velfit_dzStep = 4e-4;
    Float64 _vStep = 50.0;

    auto  result_s1 = linemodel_s1.Compute( dataStore,
                      _spc,
                      _spcContinuum,
                      tplCatalog,
                      tplCategoryList,
                      m_calibrationPath,
                      restraycatalog,
                      m_opt_linetypefilter,
                      m_opt_lineforcefilter,
                      lambdaRange,
                      _redshifts,
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
                      m_opt_firstpass_largegridstep,
                      m_opt_firstpass_largegridsampling,
                      m_opt_rigidity,
                      m_opt_tplratio_reldirpath,
                      m_opt_offsets_reldirpath,
                      m_opt_em_velocity_fit_min,
                      m_opt_em_velocity_fit_max,
                      _vStep,
                      m_opt_abs_velocity_fit_min,
                      m_opt_abs_velocity_fit_max,
                      _vStep);
    if( !result_s1 )
    {
        Log.LogInfo( "Zweimodel - Failed to compute linemodel s1");
        return false;
    }
//    else{
//        // Store linemodel chisquare results
//        dataStore.StoreScopedGlobalResult( scopeStr.c_str(), result_s1 );
//        //save linemodel fitting and spectrum-model results
//        linemodel_s1.storeGlobalModelResults(dataStore);
//        return true;
//    }

    // ---------------------------------------------------
    Log.LogInfo("Zweimodel - Computing second linemodel : S2");
    Log.LogInfo("Zweimodel - ===============================================");
    Log.LogInfo("Zweimodel - using zero continuum (probably overrided inside multimodel)");
    //WARNING: using dumb continuum zero here.
    //WARNING: Note that, in the multimodel, the continuum is re-estimated from the combined inside the operator anyway.
    CSpectrum _spcContinuum_s2 = *contSpectrum;
    _spcContinuum_s2.SetMedianWinsize(_spcContinuum.GetMedianWinsize());
    _spcContinuum_s2.SetDecompScales(_spcContinuum.GetDecompScales());
    _spcContinuum_s2.SetContinuumEstimationMethod(_spcContinuum.GetContinuumEstimationMethod());
    _spcContinuum_s2.SetWaveletsDFBinPath(_spcContinuum.GetWaveletsDFBinPath());
    CSpectrumFluxAxis spcfluxAxis2 = _spcContinuum_s2.GetFluxAxis();
    for(Int32 k=0; k<spcfluxAxis2.GetSamplesCount(); k++)
    {
        spcfluxAxis2[k] = 0.0;
    }

    COperatorLineModel linemodel_s2;
    linemodel_s2.m_secondPass_velfit_dzInfLim = 4e-4;
    linemodel_s2.m_secondPass_velfit_dzSupLim = 4e-4;
    linemodel_s2.m_secondPass_velfit_dzStep = 4e-4;
    auto  result_s2 = linemodel_s2.Compute( dataStore,
                      *contSpectrum,
                      _spcContinuum_s2,
                      tplCatalog,
                      tplCategoryList,
                      m_calibrationPath,
                      restraycatalog,
                      m_opt_linetypefilter,
                      m_opt_lineforcefilter,
                      lambdaRange,
                      _redshifts,
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
                      m_opt_firstpass_largegridstep,
                      m_opt_firstpass_largegridsampling,
                      m_opt_rigidity,
                      m_opt_tplratio_reldirpath,
                      m_opt_offsets_reldirpath,
                      m_opt_em_velocity_fit_min,
                      m_opt_em_velocity_fit_max,
                      _vStep,
                      m_opt_abs_velocity_fit_min,
                      m_opt_abs_velocity_fit_max,
                      _vStep);
    if( !result_s2 )
    {
        Log.LogInfo( "Zweimodel - Failed to compute linemodel s2");
        return false;
    }


    //get the s1 redshift candidates
    Log.LogInfo( "Zweimodel - Retrieving candidates for s1");
    std::vector<Float64> zcandidates_s1;
    auto lineModelResult_s1 = std::dynamic_pointer_cast<const CLineModelResult>( result_s1 );
    for( Int32 i=0; i<lineModelResult_s1->ExtremaResult.Extrema.size(); i++ )
    {
        Float64 zc = lineModelResult_s1->ExtremaResult.Extrema[i];
        zcandidates_s1.push_back(zc);
        Log.LogInfo( "Zweimodel - Candidate s1 - #%d : %.5f", i, zc);
    }
    //get the s2 redshift candidates
    Log.LogInfo( "Zweimodel - Retrieving candidates for s2");
    std::vector<Float64> zcandidates_s2;
    std::vector<Float64> meritcandidates_s2;
    auto lineModelResult_s2 = std::dynamic_pointer_cast<const CLineModelResult>( result_s2 );
    for( Int32 i=0; i<lineModelResult_s2->ExtremaResult.Extrema.size(); i++ )
    {
        Float64 zc = lineModelResult_s2->ExtremaResult.Extrema[i];
        zcandidates_s2.push_back(zc);
        meritcandidates_s2.push_back(lineModelResult_s2->ExtremaResult.ExtremaMerit[i]);
        Log.LogInfo( "Zweimodel - Candidate s2 - #%d : %.5f", i, zc);
    }


    std::shared_ptr<CLineModelResult> result_s1_c2;
    std::shared_ptr<CModelSpectrumResult> best_s1spcModelSpectrum;
    std::shared_ptr<CModelSpectrumResult> best_s1contModelSpectrum;
    std::shared_ptr<CModelSpectrumResult> best_s2spcModelSpectrum;
    std::shared_ptr<CModelSpectrumResult> best_s2contModelSpectrum;
    std::vector<std::vector<Float64>> combined_merits;
    Float64 bestmerit = DBL_MAX;
    Float64 bestC2z = -1.0;
    Float64 bestC1z = -1.0;
    for(Int32 kzs2=0; kzs2<zcandidates_s2.size(); kzs2++)
    {
        Log.LogInfo("Zweimodel - Computing multiroll model for s1 with contaminant from s2z%d", kzs2);
        Log.LogInfo("Zweimodel - ===============================================");
        std::vector<Float64> combined_merits_s2zX;

        COperatorLineModel linemodel_s1c2X;
        Int32 iRollContaminated_s1 = 0; //hardcoded for now
        Float64 contLambdaOffset_s2tos1 = offsetLambdaContaminant;//previously tested with hardcoded 2000A; //This is a simple modelization of the needed regrid of s2 on s1 grid (kind of due to s1/s2 radec distance)
        std::shared_ptr<CModelSpectrumResult> contModelSpectrum = linemodel_s2.GetModelSpectrumResult(kzs2);
        std::shared_ptr<CSpectraFluxResult> contModelContinuumSpectrum = linemodel_s2.GetModelSpectrumContinuumResult(kzs2);
        std::shared_ptr<CModelSpectrumResult> contModelContinuumSubtractedSpectrum = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(contModelSpectrum->GetSpectrum()) );
        for(Int32 i=0; i<contModelSpectrum->GetSpectrum().GetFluxAxis().GetSamplesCount(); i++)
        {
            contModelContinuumSubtractedSpectrum->GetSpectrum().GetFluxAxis()[i] = contModelSpectrum->GetSpectrum().GetFluxAxis()[i] - contModelContinuumSpectrum->fluxes[i];
            //contModelContinuumSubtractedSpectrum->GetSpectrum().GetFluxAxis()[i] = 0.0;
        }
        linemodel_s1c2X.initContaminant(contModelContinuumSubtractedSpectrum, iRollContaminated_s1, contLambdaOffset_s2tos1);

        Float64 _opt_twosteplargegridstep = -1;
        std::string _opt_twosteplargegridsampling = "log";
        linemodel_s1c2X.m_secondPass_extensionradius = 0.0;
        linemodel_s1c2X.m_secondPass_velfit_dzInfLim = 0; //no z refinement
        linemodel_s1c2X.m_secondPass_velfit_dzSupLim = 0; //no z refinement
        linemodel_s1c2X.m_secondPass_velfit_dzStep = 2e-4;
        Int32 _opt_extremacount = -1; //no extrema search calculating on all redshifts
        auto result_s1_c2zX = linemodel_s1c2X.Compute( dataStore,
                                                     _spc,
                                                     _spcContinuum,
                                                     tplCatalog,
                                                     tplCategoryList,
                                                     m_calibrationPath,
                                                     restraycatalog,
                                                     m_opt_linetypefilter,
                                                     m_opt_lineforcefilter,
                                                     lambdaRange,
                                                     zcandidates_s1,
                                                     _opt_extremacount,
                                                     m_opt_fittingmethod,
                                                     m_opt_continuumcomponent,
                                                     m_opt_lineWidthType,
                                                     m_opt_resolution,
                                                     m_opt_velocity_emission,
                                                     m_opt_velocity_absorption,
                                                     m_opt_continuumreest,
                                                     m_opt_rules,
                                                     m_opt_velocityfit,
                                                     _opt_twosteplargegridstep,
                                                     _opt_twosteplargegridsampling,
                                                     m_opt_rigidity,
                                                     m_opt_tplratio_reldirpath,
                                                     m_opt_offsets_reldirpath,
                                                     m_opt_em_velocity_fit_min,
                                                     m_opt_em_velocity_fit_max,
                                                     _vStep,
                                                     m_opt_abs_velocity_fit_min,
                                                     m_opt_abs_velocity_fit_max,
                                                     _vStep);

        if( !result_s1_c2zX )
        {
            Log.LogInfo( "Zweimodel - Failed to compute linemodel s1_c2z%d", kzs2);
            return false;
        }
        else{
            //now computing s2zX with c1zY
            auto lineModelResult_s1c2zX = std::dynamic_pointer_cast<const CLineModelResult>( result_s1_c2zX );

            for(Int32 kzs1=0; kzs1<zcandidates_s1.size(); kzs1++)
            {

                Log.LogInfo("Zweimodel - Computing multiroll model for s2 with contaminant from c1zY%d", kzs1);
                Log.LogInfo("Zweimodel - ===============================================");

                Int32 idx_linemodelExtrema_s1c2zX = -1;
                for(Int32 i=0; i<lineModelResult_s1c2zX->ExtremaResult.Extrema.size(); i++)
                {
                    if(zcandidates_s1[kzs1]==lineModelResult_s1c2zX->ExtremaResult.Extrema[i])
                    {
                        idx_linemodelExtrema_s1c2zX = i;
                    }
                }
                Float64 redshift1 = lineModelResult_s1c2zX->ExtremaResult.Extrema[idx_linemodelExtrema_s1c2zX];
                Float64 merits1 = lineModelResult_s1c2zX->ExtremaResult.ExtremaMerit[idx_linemodelExtrema_s1c2zX];

                COperatorLineModel linemodel_s2c1Y;
                Int32 iRollContaminated_s2 = 0; //hardcoded for now
                Float64 contLambdaOffset_s1tos2 = -offsetLambdaContaminant;
                std::shared_ptr<CModelSpectrumResult> contModelSpectrum = linemodel_s1c2X.GetModelSpectrumResult(idx_linemodelExtrema_s1c2zX);
                std::shared_ptr<CSpectraFluxResult> contModelContinuumSpectrum = linemodel_s1c2X.GetModelSpectrumContinuumResult(idx_linemodelExtrema_s1c2zX);
                std::shared_ptr<CModelSpectrumResult> contModelContinuumSubtractedSpectrum = std::shared_ptr<CModelSpectrumResult>( new CModelSpectrumResult(contModelSpectrum->GetSpectrum()) );
                for(Int32 i=0; i<contModelSpectrum->GetSpectrum().GetFluxAxis().GetSamplesCount(); i++)
                {
                    contModelContinuumSubtractedSpectrum->GetSpectrum().GetFluxAxis()[i] = contModelSpectrum->GetSpectrum().GetFluxAxis()[i] - contModelContinuumSpectrum->fluxes[i];
                    //contModelContinuumSubtractedSpectrum->GetSpectrum().GetFluxAxis()[i] = 0.0;
                }
                linemodel_s2c1Y.initContaminant(contModelContinuumSubtractedSpectrum, iRollContaminated_s2, contLambdaOffset_s1tos2);

                Float64 _opt_twosteplargegridstep = -1;
                std::string _opt_twosteplargegridsampling = "log";
                linemodel_s2c1Y.m_secondPass_extensionradius = 0.0;
                linemodel_s2c1Y.m_secondPass_velfit_dzInfLim = 0;
                linemodel_s2c1Y.m_secondPass_velfit_dzSupLim = 0;
                linemodel_s2c1Y.m_secondPass_velfit_dzStep = 2e-4;
                Int32 _opt_extremacount = -1; //no extrema search calculating on all redshifts
                std::vector<Float64> redshifts_s2(1, zcandidates_s2[kzs2]);
                auto result_s2_c1zY = linemodel_s2c1Y.Compute( dataStore,
                                                            *contSpectrum,
                                                            _spcContinuum_s2,
                                                             tplCatalog,
                                                             tplCategoryList,
                                                             m_calibrationPath,
                                                             restraycatalog,
                                                             m_opt_linetypefilter,
                                                             m_opt_lineforcefilter,
                                                             lambdaRange,
                                                             redshifts_s2,
                                                             _opt_extremacount,
                                                             m_opt_fittingmethod,
                                                             m_opt_continuumcomponent,
                                                             m_opt_lineWidthType,
                                                             m_opt_resolution,
                                                             m_opt_velocity_emission,
                                                             m_opt_velocity_absorption,
                                                             m_opt_continuumreest,
                                                             m_opt_rules,
                                                             m_opt_velocityfit,
                                                             _opt_twosteplargegridstep,
                                                             _opt_twosteplargegridsampling,
                                                             m_opt_rigidity,
                                                             m_opt_tplratio_reldirpath,
                                                             m_opt_offsets_reldirpath,
                                                             m_opt_em_velocity_fit_min,
                                                             m_opt_em_velocity_fit_max,
                                                             _vStep,
                                                             m_opt_abs_velocity_fit_min,
                                                             m_opt_abs_velocity_fit_max,
                                                             _vStep);

                if( !result_s2_c1zY )
                {
                    Log.LogInfo( "Zweimodel - Failed to compute linemodel s2_c1z%d", kzs1);
                    return false;
                }else{

                    auto lineModelResult_s2c1zY = std::dynamic_pointer_cast<const CLineModelResult>( result_s2_c1zY );
                    Float64 redshift2 = lineModelResult_s2c1zY->ExtremaResult.Extrema[0];
                    Float64 merits2 = lineModelResult_s2c1zY->ExtremaResult.ExtremaMerit[0];
                    Float64 combinedMerit = merits2 + merits1;
                    combined_merits_s2zX.push_back(combinedMerit);

                    if(bestmerit>combinedMerit)
                    {
                        bestC2z = redshift2;
                        bestC1z = redshift1;
                        bestmerit = combinedMerit;
                        result_s1_c2 = std::dynamic_pointer_cast<CLineModelResult>( result_s1_c2zX );
                        best_s1spcModelSpectrum = linemodel_s1c2X.GetModelSpectrumResult(idx_linemodelExtrema_s1c2zX);
                        best_s1contModelSpectrum = linemodel_s1c2X.GetContaminantSpectrumResult();
                        best_s2spcModelSpectrum = linemodel_s2c1Y.GetModelSpectrumResult(0);
                        best_s2contModelSpectrum = linemodel_s2c1Y.GetContaminantSpectrumResult();
                    }
                }



            }
            combined_merits.push_back(combined_merits_s2zX);
        }
    }

    if( !result_s1_c2 )
    {
        Log.LogInfo( "Zweimodel - Failed to compute linemodel s1_c2");
        return false;
    }
    else{
        Log.LogInfo( "Zweimodel - Successfully computed linemodel s1_c2: best S1_z=%.5f, best S2_z=%.5f", bestC1z, bestC2z);
        // Store linemodel chisquare results
        dataStore.StoreScopedGlobalResult( scopeStr.c_str(), result_s1_c2 );
        //save linemodel fitting and spectrum-model results
        std::string fname_s1spc = (boost::format("linemodel_s1spc_extrema_0")).str();
        dataStore.StoreScopedGlobalResult( fname_s1spc.c_str(), best_s1spcModelSpectrum );
        std::string fname_s1cont = (boost::format("linemodel_s1cont_extrema_0")).str();
        dataStore.StoreScopedGlobalResult( fname_s1cont.c_str(), best_s1contModelSpectrum );
        std::string fname_s2spc = (boost::format("linemodel_s2spc_extrema_0")).str();
        dataStore.StoreScopedGlobalResult( fname_s2spc.c_str(), best_s2spcModelSpectrum );
        std::string fname_s2cont = (boost::format("linemodel_s2cont_extrema_0")).str();
        dataStore.StoreScopedGlobalResult( fname_s2cont.c_str(), best_s2contModelSpectrum );

        std::shared_ptr<CZweiModelResult> zweimodelResult = std::shared_ptr<CZweiModelResult>( new CZweiModelResult(zcandidates_s1, zcandidates_s2, combined_merits));
        std::string fname_result = (boost::format("zweimodel")).str();
        dataStore.StoreScopedGlobalResult( fname_result.c_str(), zweimodelResult );
    }

    return true;
}

Int32 CZweiModelSolve::getContaminantNameFromFile( const char* filePath, std::string spcid, Int32 colID,
                                                   std::string& contSpcFileName,
                                                   std::string& contErrorFileName,
                                                   Float64& offsetLambdaContaminant,
                                                   Int32 reverseInclusion )
{    
    contSpcFileName = "";
    contErrorFileName = "";

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
            Int32 nskip = colID;
            for(Int32 i=0; i<nskip; i++)
            {
                ++it;
            }
            if( it != tok.end() )
            {

                try
                {
                    contSpcFileName = *it;

                    size_t pos = contSpcFileName.find( "_F" );
                    if ( pos != string::npos ) {
                        contErrorFileName = contSpcFileName;
                        contErrorFileName.replace( pos, 2, "_ErrF" );
                    }else{
                        contErrorFileName = "";
                        Log.LogError( "Failed to find spc tag in contaminant name" );
                        return false;
                    }

                    //now read the lambda offset
                    ++it;
                    if( it != tok.end() )
                    {
                        try
                        {
                            offsetLambdaContaminant = lexical_cast<double>(*it);
                        }
                        catch (bad_lexical_cast&)
                        {
                            return false;
                        }
                    }else{
                        Log.LogError( "Failed to find offsetLambdaContaminant value for the given contaminant spectrum" );
                        return false;
                    }

                    file.close();
                    return true;
                }
                catch (bad_lexical_cast&)
                {
                    return false;
                }
            }

        }
    }
    file.close();

    return true;
}

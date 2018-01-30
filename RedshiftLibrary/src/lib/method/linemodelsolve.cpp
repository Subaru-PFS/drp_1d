#include <RedshiftLibrary/method/linemodelsolve.h>

#include <RedshiftLibrary/log/log.h>

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
    desc.append("\tparam: linemodel.rigidity = {""rules"", ""tplcorr"", ""tplshape""}\n");
    desc.append("\tparam: linemodel.linewidthtype = {""instrumentdriven"", ""velocitydriven"",  ""combined"",  ""nispvsspsf201707"", ""fixed""}\n");
    desc.append("\tparam: linemodel.instrumentresolution = <float value>\n");
    desc.append("\tparam: linemodel.velocityemission = <float value>\n");
    desc.append("\tparam: linemodel.velocityabsorption = <float value>\n");
    desc.append("\tparam: linemodel.continuumreestimation = {""no"", ""onlyextrema"", ""always""}\n");
    desc.append("\tparam: linemodel.rules = {""all"", ""balmer"", ""strongweak"", ""superstrong"", ""ratiorange"", ""ciiiratio"", ""no""}\n");
    desc.append("\tparam: linemodel.extremacount = <float value>\n");
    desc.append("\tparam: linemodel.velocityfit = {""yes"", ""no""}\n");
    desc.append("\tparam: linemodel.emvelocityfitmin = <float value>\n");
    desc.append("\tparam: linemodel.emvelocityfitmax = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitmin = <float value>\n");
    desc.append("\tparam: linemodel.absvelocityfitmax = <float value>\n");
    desc.append("\tparam: linemodel.fastfitlargegridstep = <float value>, deactivated if negative or zero\n");
    desc.append("\tparam: linemodel.pdfcombination = {""marg"", ""bestchi2""}\n");
    desc.append("\tparam: linemodel.stronglinesprior = {""yes"", ""no""}\n");
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
    dataStore.GetScopedParam( "linemodel.fastfitlargegridstep", m_opt_twosteplargegridstep, 0.001 );
    dataStore.GetScopedParam( "linemodel.continuumcomponent", m_opt_continuumcomponent, "fromspectrum" );
    dataStore.GetScopedParam( "linemodel.rigidity", m_opt_rigidity, "rules" );
    dataStore.GetScopedParam( "linemodel.linewidthtype", m_opt_lineWidthType, "velocitydriven" );
    dataStore.GetScopedParam( "linemodel.instrumentresolution", m_opt_resolution, 2350.0 );
    dataStore.GetScopedParam( "linemodel.velocityemission", m_opt_velocity_emission, 100.0 );
    dataStore.GetScopedParam( "linemodel.velocityabsorption", m_opt_velocity_absorption, 300.0 );
    dataStore.GetScopedParam( "linemodel.velocityfit", m_opt_velocityfit, "yes" );
    if(m_opt_velocityfit=="yes"){
        dataStore.GetScopedParam( "linemodel.emvelocityfitmin", m_opt_em_velocity_fit_min, 20.0 );
        dataStore.GetScopedParam( "linemodel.emvelocityfitmax", m_opt_em_velocity_fit_max, 500.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitmin", m_opt_abs_velocity_fit_min, 150.0 );
        dataStore.GetScopedParam( "linemodel.absvelocityfitmax", m_opt_abs_velocity_fit_max, 500.0 );
    }
    dataStore.GetScopedParam( "linemodel.continuumreestimation", m_opt_continuumreest, "no" );
    dataStore.GetScopedParam( "linemodel.rules", m_opt_rules, "all" );
    dataStore.GetScopedParam( "linemodel.extremacount", m_opt_extremacount, 10.0 );
    dataStore.GetScopedParam( "linemodel.stronglinesprior", m_opt_stronglinesprior, "no");
    dataStore.GetScopedParam( "linemodel.pdfcombination", m_opt_pdfcombination, "marg");
    dataStore.GetScopedParam( "linemodel.saveintermediateresults", m_opt_saveintermediateresults, "no");

    //Auto-correct fitting method
    //std::string forcefittingmethod = "ones";
    std::string forcefittingmethod = "individual";
    if(m_opt_rigidity=="tplshape" && m_opt_fittingmethod != forcefittingmethod)
    {
        m_opt_fittingmethod = forcefittingmethod;
        dataStore.SetScopedParam("linemodel.fittingmethod", m_opt_fittingmethod);
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
    }

    Log.LogInfo( "    -continuumcomponent: %s", m_opt_continuumcomponent.c_str());
    Log.LogInfo( "    -continuumreestimation: %s", m_opt_continuumreest.c_str());
    Log.LogInfo( "    -extremacount: %.3f", m_opt_extremacount);
    Log.LogInfo( "    -fastfitlargegridstep: %.6f", m_opt_twosteplargegridstep);

    Log.LogInfo( "    -pdf-stronglinesprior: %s", m_opt_stronglinesprior.c_str());
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
								       const TFloat64List& redshifts )
{
    CDataStore::CAutoScope resultScope( dataStore, "linemodelsolve" );

    PopulateParameters( dataStore );
    Int32 retSolve = Solve( dataStore, spc, spcWithoutCont, tplCatalog, tplCategoryList, restraycatalog, lambdaRange, redshifts);

    if(retSolve){

        /* ------------------------  COMPUTE POSTMARG PDF  --------------------------  */
        Log.LogInfo("linemodelsolve: Pdfz computation");

        std::string scope = "linemodelsolve.linemodel";
        auto results = dataStore.GetGlobalResult( scope.c_str() );
        if(results.expired())
        {
            Log.LogError("linemodelsolve: Unable to retriev linemodel results");
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
                result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0);
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
                result_chisquaretplshape->Init( result->Redshifts, result->restRayList, 0);
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


Int32 CLineModelSolve::CombinePDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result, std::string opt_rigidity, std::string opt_combine, std::string opt_stronglinesprior)
{
    //hardcoded prior parameter
    bool zPriorStrongLinePresence = (opt_stronglinesprior=="yes");
    if(zPriorStrongLinePresence)
    {
        Log.LogInfo("Linemodel: Pdfz computation: StrongLinePresence prior enabled");
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
            zPrior->valProbaLog = pdfz.GetStrongLinePresenceLogZPrior(strongLinePresence);
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
        std::vector<TFloat64List> priorsTplshapes;
        for(Int32 k=0; k<result->ChiSquareTplshapes.size(); k++)
        {
            TFloat64List _prior;
            if(zPriorStrongLinePresence)
            {
                std::vector<bool> strongLinePresence = result->StrongELPresentTplshapes[k];
                _prior = pdfz.GetStrongLinePresenceLogZPrior(strongLinePresence);
            }else
            {
                _prior = pdfz.GetConstantLogZPrior(result->Redshifts.size());
            }
            priorsTplshapes.push_back(_prior);
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
            retPdfz = pdfz.Marginalize( result->Redshifts, ChiSquareTplshapesCorrected, priorsTplshapes, cstLog, postmargZResult);
        }else{
            retPdfz = pdfz.BestProba( result->Redshifts, ChiSquareTplshapesCorrected, priorsTplshapes, cstLog, postmargZResult);
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
                catch (bad_lexical_cast)
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
                catch (bad_lexical_cast)
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
Int32 getValueFromRefFile( const char* filePath, std::string spcid, Int32 colID, Float64& zref, Int32 reverseInclusion )
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
                catch (bad_lexical_cast)
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

    //Hack: load the zref values
    std::vector<Float64> _redshifts;
    if(false)
    {
        Log.LogInfo( "Linemodel - hacking zref enabled");
        Float64 zref = -1.0;
        namespace fs = boost::filesystem;
        Int32 reverseInclusionForIdMatching = 0;
        //*
        // Euclid case for process at zref
        fs::path refFilePath("/home/aschmitt/amazed_cluster/datasets/euclid/euclidsim2016/EUC-TEST-TUGALSPC-2016-03_export20170302_3z4mag4lfhabins_TrueFlux/reference_correctedFastSim.list");
        //fs::path refFilePath("/sps/euclid/Users/schmitt/amazed_cluster/amazed_cluster/datasets/euclid/euclidsim2016/EUC-TEST-TUGALSPC-2016-03_export20170302_3z4mag4lfhabins_TrueFlux/reference_correctedFastSim.list");
        Int32 substring_start = 0;
        Int32 substring_n = spc.GetName().size();
        Int32 colId = 2; //starts at 1, so that id_column is usually 1
        reverseInclusionForIdMatching = 1;
        //*/
        /*
        // SDSS case for process at zref
        //fs::path refFilePath("/home/aschmitt/amazed_cluster/datasets/sdss/sdss_201707/reference_SDSS_spectra_bg10k.txt");
        //fs::path refFilePath("/sps/euclid/Users/schmitt/amazed_cluster/datasets/sdss/sdss_201707/reference_SDSS_spectra_bg10k.txt");
        fs::path refFilePath("/home/aschmitt/data/sdss/sdss_201707/reference_SDSS_spectra_bg10k.txt");
        Int32 substring_start = 5;
        Int32 substring_n = 15;
        Int32 colId = 2; //starts at 1, so that id_column is usually 1
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

        if(true) //computing only on zref, or on a zrange around zref
        {
            Float64 deltaZrangeHalf = 0.5e-2; //override zrange
            Float64 stepZ = 1e-4;
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
        m_opt_twosteplargegridstep = 0; //override fastfitlargegrid
        Log.LogInfo( "Linemodel - hack - Loaded zref for spc %s : zref=%f", spc.GetName().c_str(), zref);
    }else{
        _redshifts = redshifts;
    }

    // Compute merit function
    COperatorLineModel linemodel;
    auto  result = linemodel.Compute( dataStore,
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
                      m_opt_twosteplargegridstep,
                      m_opt_rigidity,
                      m_opt_em_velocity_fit_min,
                      m_opt_em_velocity_fit_max,
                      m_opt_abs_velocity_fit_min,
                      m_opt_abs_velocity_fit_max);

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


#include <RedshiftLibrary/processflow/processflow.h>
#include <RedshiftLibrary/processflow/autoscope.h>

#include <RedshiftLibrary/linemodel/calibrationconfig.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>

#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/method/templatefittingsolve.h>
#include <RedshiftLibrary/method/linematchingsolve.h>
#include <RedshiftLibrary/method/tplcombinationsolve.h>
#include <RedshiftLibrary/method/linemodelsolve.h>
#include <RedshiftLibrary/method/classificationsolve.h>
#include <RedshiftLibrary/method/classificationresult.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>
#include <RedshiftLibrary/reliability/zqual.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <cstdio>
#include <cfloat>

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

    TFloat64Range lambdaRange;
    TFloat64Range redshiftRange;
    Float64       redshiftStep;
    Float64       maxCount; 
    Float64       redshiftseparation;

    ctx.GetParameterStore()->Get( "lambdarange", lambdaRange );
    ctx.GetParameterStore()->Get( "redshiftrange", redshiftRange );
    ctx.GetParameterStore()->Get( "redshiftstep", redshiftStep );
    ctx.GetParameterStore()->Get( "extremaredshiftseparation", redshiftseparation, 2*0.005);//todo: decide on the default values using latest analyses plots



    //Override z-search grid for line measurement: load the zref values from a tsv catalog file (col0: spc name, col1: zref float value)
    // std::string opt_linemeas_catalog_path;
    //ctx.GetParameterStore()->Get( "linemeascatalog", opt_linemeas_catalog_path, "" );
    // TFloat64List redshifts=raw_redshifts;
    //Float64 zref = -1.0;
    /*
    if(opt_linemeas_catalog_path!="")
    {
        Log.LogInfo( "Override z-search range !");
        Int32 reverseInclusionForIdMatching = 0; //0: because the names must match exactly, but: linemeas catalog includes the extension (.fits) and spc.GetName doesn't.
    
        Int32 colId = 2;//starts at 1, so that (for the linemeas_catalog) id_column=1, zref_column=2
        bfs::path refFilePath(opt_linemeas_catalog_path.c_str());
        
        if ( bfs::exists(refFilePath) )
          {
            std::string spcSubStringId = spectrumName;
            Log.LogInfo( "Override z-search: using spc-string: %s", spcSubStringId.c_str());
            getValueFromRefFile( refFilePath.c_str(), spcSubStringId, zref, reverseInclusionForIdMatching);
        }
        
        if(zref==-1.0)
        {
          throw std::runtime_error("Override z-search: unable to find zref!");
    }


        Float64 stepZ = 1e-5;
        Float64 deltaZrangeHalf = -1;//0.5e-2; //override zrange
        ctx.GetParameterStore()->Get( "linemeas.dzhalf", deltaZrangeHalf, -1);
        ctx.GetParameterStore()->Get( "linemeas.dzstep", stepZ, 1e-5);
        bool computeOnZrange=false;
        if(deltaZrangeHalf>0.0)
    {
            computeOnZrange=true;
        }
        if(computeOnZrange) //computing only on zref, or on a zrange around zref
        {
            //make sure to include zref in redshift vector
            Float64 nStepsZhalf = std::ceil(deltaZrangeHalf/stepZ);
            redshifts.resize(nStepsZhalf*2 + 1); 
            for(Int32 kz= 0; kz<nStepsZhalf +1 ; kz++)
            {
                Float64 _z = zref + kz*stepZ;
                redshifts[nStepsZhalf + kz] = _z;
        }
            for(Int32 kz= 1; kz<nStepsZhalf + 1; kz++)
        {
                Float64 _z = zref - kz*stepZ;
                redshifts[nStepsZhalf - kz]= _z;
        }
            Log.LogInfo( "Override z-search: zmin=%.5f, zmax=%.5f, zstep=%.5f", redshifts[0], redshifts[redshifts.size()-1], stepZ);
        }else{
            redshifts.push_back(zref);
        }

        if(methodName=="linemodel"){
          
            ctx.GetDataStore().SetScopedParam( "linemodelsolve.linemodel.extremacount", 1.0);
            Log.LogInfo( "Override z-search: Using overriden linemodelsolve.linemodel.extremacount: %f", 1.0);
            ctx.GetDataStore().SetScopedParam( "linemodelsolve.linemodel.firstpass.largegridstep", stepZ);

            Log.LogInfo( "Override z-search: Using overriden linemodelsolve.linemodel.firstpass.largegridstep: %f", stepZ);
            Log.LogInfo( "Override z-search: Using overriden half zrange around zref: %f", deltaZrangeHalf);
            Log.LogInfo( "Override z-search: Using overriden dzstep: %f", stepZ);
        }

        Log.LogInfo( "Override z-search: Using overriden zref for spc %s : zref=%f", spectrumName.c_str(), zref);
    }else{
        redshifts = raw_redshifts;
    }
    ctx.GetDataStore().SetScopedParam( "linemodelsolve.linemodel.zref", zref);

    if(redshifts.size() < 1)
    {
      Log.LogError("Unable to initialize the z-search grid (size=%d)", redshifts.size());
      throw std::runtime_error("Unable to initialize the z-search grid");
    }


    */
    //retrieve the calibration dir path
    std::string calibrationDirPath = ctx.GetParameterStore()->Get<std::string>( "calibrationDir");

    //************************************
    // Stellar method


    if(ctx.GetParameterStore()->Get<std::string>( "enablestellarsolve")=="yes"){
      //      CAutoScope resultScope( ctx.m_ScopeStack, "stellarsolve" );

        Log.LogInfo("Processing stellar fitting");
        CMethodTemplateFittingSolve solve(ctx.m_ScopeStack,"star");
        solve.Compute(ctx.GetInputContext(),
                      ctx.GetResultStore(),
                      ctx.m_ScopeStack);
   
    }

    // Quasar method
    //    std::shared_ptr<CSolveResult> qsoResult;
    std::string enableQsoFitting;
    ctx.GetParameterStore()->Get( "enableqsosolve", enableQsoFitting, "no" );
    Log.LogInfo( "QSO solve enabled : %s", enableQsoFitting.c_str());
    /*
    if(enableQsoFitting=="yes"){
        CDataStore::CAutoScope resultScope( ctx.m_ScopeStack, "qsosolve" );

        std::string qso_method;
        ctx.GetParameterStore()->Get( "qsosolve.method", qso_method, "templatefittingsolve" );
        Float64 overlapThreshold;
        ctx.GetParameterStore()->Get( "qsosolve.overlapThreshold", overlapThreshold, 1.0);
        std::string opt_spcComponent;
        ctx.GetParameterStore()->GetScopedParam( "qsosolve.spectrum.component", opt_spcComponent, "raw" );
        std::string opt_interp;
        ctx.GetParameterStore()->GetScopedParam( "qsosolve.interpolation", opt_interp, "precomputedfinegrid" );
        std::string opt_extinction;
        ctx.GetParameterStore()->GetScopedParam( "qsosolve.extinction", opt_extinction, "yes" );
        std::string opt_dustFit;
        ctx.GetParameterStore()->GetScopedParam( "qsosolve.dustfit", opt_dustFit, "no" );

        // prepare the unused masks
        std::vector<CMask> maskList;
        //define the redshift search grid
        TFloat64Range qsoRedshiftRange=TFloat64Range(0.0, 6.0);
        Float64 qsoRedshiftStep = 5e-4;
        Log.LogInfo("QSO fitting redshift range = [%.5f, %.5f], step=%.6f", qsoRedshiftRange.GetBegin(), qsoRedshiftRange.GetEnd(), qsoRedshiftStep);
        TFloat64List qso_redshifts;
        if(redshiftSampling=="log")
        {
            qso_redshifts = qsoRedshiftRange.SpreadOverLog( qsoRedshiftStep );
        }else{
            qso_redshifts = qsoRedshiftRange.SpreadOver( qsoRedshiftStep );
        }
        DebugAssert( qso_redshifts.size() > 0 );

        Log.LogInfo("Processing QSO fitting");
        if(qso_method=="templatefittingsolve"){
            CMethodTemplateFittingSolve solve;
            qsoResult = solve.Compute( ctx.GetDataStore(),
                                       ctx.GetSpectrum(),
                                       ctx.GetTemplateCatalog(),
                                       ctx.GetQSOCategoryList(),
                                       ctx.GetInputContext().m_lambdaRange,
                                       qso_redshifts,
                                       overlapThreshold,
                                       maskList,
                                       "qso_zPDF",
                                       redshiftseparation,
                                       opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);
        
        }else if(qso_method=="tplcombinationsolve"){
            opt_interp="lin";
            opt_extinction="no";
            opt_dustFit="no";
            CMethodTplcombinationSolve solve;
            qsoResult = solve.Compute( ctx.GetDataStore(),
                                       ctx.GetSpectrum(),
                                       ctx.GetTemplateCatalog(),
                                       ctx.GetQSOCategoryList(),
                                       ctx.GetInputContext().m_lambdaRange,
                                       qso_redshifts,
                                       overlapThreshold,
                                       maskList,
                                       "qso_zPDF",
                                       redshiftseparation,
                                       opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);
        }else if(qso_method=="linemodel"){
            Log.LogInfo("Linemodel qso fitting...");
            CLineModelSolve Solve(calibrationDirPath);
            qsoResult = Solve.Compute( ctx.GetDataStore(),
                                       ctx.GetSpectrum(),
                                       ctx.GetTemplateCatalog(),
                                       ctx.GetQSOCategoryList(),
                                       ctx.GetRayCatalog(),
                                       ctx.GetInputContext().m_lambdaRange,
                                       qso_redshifts,
                                       "qso_zPDF",
                                       redshiftseparation);
        }else{
            throw std::runtime_error("Problem found while parsing the qso method parameter");
        }

        //finally save the qso fitting results
        if( qsoResult ) {
            Log.LogInfo("Saving qso fitting results");
            ctx.GetDataStore().StoreScopedGlobalResult( "qsoresult", qsoResult );
        }else{
            Log.LogError( "Unable to store qso result.");
        }
        //qsoResult->preSave(ctx.GetDataStore());
    }
*/
    // Galaxy method


    if(true)
      {
        std::string methodName;
        ctx.GetParameterStore()->Get( "galaxy.method", methodName );
        boost::algorithm::to_lower(methodName);

        std::string galaxy_method_pdf_reldir = "zPDF";
        if(methodName  == "linemodelsolve" ){

          CLineModelSolve Solve(ctx.m_ScopeStack,"galaxy",calibrationDirPath);
          Solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);
          /*
            }else if(methodName  == "zweimodelsolve" ){

            CZweiModelSolve Solve(calibrationDirPath);
            mResult = Solve.Compute( ctx.GetDataStore(),
            ctx.GetSpectrum(),
            ctx.GetTemplateCatalog(),
            ctx.GetGalaxyCategoryList(),
            ctx.GetRayCatalog(),
            ctx.GetInputContext().m_lambdaRange,
            redshifts );
          */
        }else if(methodName  == "templatefittingsolve" ){
          CMethodTemplateFittingSolve solve(ctx.m_ScopeStack,"galaxy");
          solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);

        }else if(methodName  == "tplcombinationsolve" ){
          Float64 overlapThreshold;
          ctx.GetParameterStore()->Get( "tplcombinationsolve.overlapThreshold", overlapThreshold, 1.0);
          std::string opt_spcComponent;
          ctx.GetParameterStore()->GetScopedParam( "tplcombinationsolve.spectrum.component", opt_spcComponent, "raw" );
          std::string opt_interp="lin";
          ctx.GetParameterStore()->GetScopedParam( "tplcombinationsolve.interpolation", opt_interp, "lin" );
          std::string opt_extinction="no";
          //ctx.GetParameterStore()->GetScopedParam( "tplcombinationsolve.extinction", opt_extinction, "no" );
          std::string opt_dustFit="no";
          //ctx.GetDataStore().GetScopedParam( "tplcombinationsolve.dustfit", opt_dustFit, "no" );

          // prepare the unused masks
          std::vector<CMask> maskList;
          /*CMethodTplcombinationSolve solve;
            mResult = solve.Compute( ctx.GetDataStore(),
            ctx.GetSpectrum(),
            ctx.GetTemplateCatalog(),
            ctx.GetGalaxyCategoryList(),
            ctx.GetInputContext().m_lambdaRange,
            redshifts,
            overlapThreshold,
            maskList,
            galaxy_method_pdf_reldir,
            redshiftseparation,
            opt_spcComponent, opt_interp, opt_extinction, opt_dustFit);
          */

        }
        /*
          }else if(methodName  == "linematching" ){
          COperatorLineMatchingSolve Solve;
          mResult = Solve.Compute(ctx.GetDataStore(),
          ctx.GetSpectrum(),
          lambdaRange,
          redshiftRange,
          redshiftStep,
          ctx.GetRayCatalog() );

        */
        /*
          else if(methodName  == "reliability" ){
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

          }*/
        else{
          throw std::runtime_error("Problem found while parsing the method parameter");
        }
      }
    
    //mResult->preSave(ctx.GetDataStore());
    

    //estimate star/galaxy/qso classification
    Log.LogInfo("===============================================");

    CClassificationSolve classifier(ctx.m_ScopeStack,"classification");
    classifier.Compute(ctx.GetInputContext(),
                       ctx.GetResultStore(),
                       ctx.m_ScopeStack);



}


/**
 * @brief isPdfValid
 * @return
 */
Bool CProcessFlow::isPdfValid(CProcessFlowContext& ctx) const
{
    std::string scope_res = "zPDF/logposterior.logMargP_Z_data";
    auto results_pdf =  ctx.GetResultStore()->GetGlobalResult( scope_res.c_str() );
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
Int32 CProcessFlow::getValueFromRefFile( const char* filePath, std::string spcid, Float64& zref, Int32 reverseInclusion )
{
    int colID = -1;

    ifstream file;

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
        ttokenizer::iterator it;

        if (colID == -1) {
            int count = 0;
            for ( it = tok.begin(); it != tok.end() ; ++count, ++it );
            switch (count) {
            case 2:
                // Two-columns style redshift reference catalog
                colID = 2;
                break;
            case 13:
                // redshift.csv style redshift reference catalog
                colID = 3;
                break;
            default:
                Log.LogError("Invalid number of columns in reference catalog (%d)", count);
                throw std::runtime_error("Invalid number of columns in reference catalog");
            }
        }
        it = tok.begin();

        // Check if it's not a comment
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

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
#include <RedshiftLibrary/method/reliabilitysolve.h>
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

    Float64       maxCount; 
    Float64       redshiftseparation;

    ctx.GetParameterStore()->Get( "extremaredshiftseparation", redshiftseparation, 2*0.005);//todo: decide on the default values using latest analyses plots

    //retrieve the calibration dir path
    std::string calibrationDirPath = ctx.GetParameterStore()->Get<std::string>( "calibrationDir");

    //************************************
    // Stellar method


    if(ctx.GetParameterStore()->Get<std::string>( "enablestellarsolve")=="yes"){
        Log.LogInfo("Processing stellar fitting");
        CMethodTemplateFittingSolve solve(ctx.m_ScopeStack,"star");
        solve.Compute(ctx.GetInputContext(),
                      ctx.GetResultStore(),
                      ctx.m_ScopeStack);
   
    }

    // Quasar method
    //    std::shared_ptr<CSolveResult> qsoResult;
    std::string enableQsoFitting = ctx.GetParameterStore()->Get<std::string>( "enableqsosolve");
    Log.LogInfo( "QSO solve enabled : %s", enableQsoFitting.c_str());

    if(enableQsoFitting=="yes"){
        Log.LogInfo("Processing QSO fitting");
        std::string qso_method = ctx.GetParameterStore()->Get<std::string>( "qso.method");
        boost::algorithm::to_lower(qso_method);
        if(qso_method=="templatefittingsolve"){
          CMethodTemplateFittingSolve solve(ctx.m_ScopeStack,"qso");
          solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);
        }else if(qso_method=="tplcombinationsolve"){
            CMethodTplcombinationSolve solve(ctx.m_ScopeStack,"qso");
            solve.Compute(ctx.GetInputContext(),
                          ctx.GetResultStore(),
                          ctx.m_ScopeStack);
                                      }
        
        else if(qso_method=="linemodel"){
            Log.LogInfo("Linemodel qso fitting...");
            CLineModelSolve Solve(ctx.m_ScopeStack,"qso",calibrationDirPath);
            Solve.Compute( ctx.GetInputContext(),
                           ctx.GetResultStore(),
                           ctx.m_ScopeStack);
        }else{
            throw std::runtime_error("Problem found while parsing the qso method parameter");
        }
    }

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
         CMethodTplcombinationSolve solve(ctx.m_ScopeStack,"galaxy");
          solve.Compute( ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);
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
        else{
          throw std::runtime_error("Problem found while parsing the method parameter");
        }
      }
    
    
    //estimate star/galaxy/qso classification
    Log.LogInfo("===============================================");

    {
      CClassificationSolve classifier(ctx.m_ScopeStack,"classification");
      classifier.Compute(ctx.GetInputContext(),
                         ctx.GetResultStore(),
                         ctx.m_ScopeStack);
    }
    {
      CReliabilitySolve reliability(ctx.m_ScopeStack,"reliability");
      reliability.Compute(ctx.GetInputContext(),
                          ctx.GetResultStore(),
                          ctx.m_ScopeStack);
    }
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

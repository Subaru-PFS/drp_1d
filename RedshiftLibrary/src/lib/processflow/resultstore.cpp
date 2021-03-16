#include <RedshiftLibrary/processflow/resultstore.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/statistics/pdfcandidateszresult.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/common/exception.h>
#include <RedshiftLibrary/common/formatter.h>

using namespace NSEpic;

namespace bfs = boost::filesystem;


COperatorResultStore::COperatorResultStore(const TScopeStack& scope):
  CScopeStore(scope)
{

}

COperatorResultStore::~COperatorResultStore()
{

}

void COperatorResultStore::StoreResult( TResultsMap& map, const std::string& path, const std::string& name,
                                        std::shared_ptr<const COperatorResult> result )
{
    std::string scopedName;
    if( ! path.empty() ) {
        scopedName = path;
        scopedName.append( "." );
    }
    scopedName.append( name );

    TResultsMap::iterator it = map.find( name );
    if( it != map.end() )
    {
        DebugError( "Result already exist");
        return;
    }

    map[ scopedName ] = result;
}


void COperatorResultStore::StorePerTemplateResult( const CTemplate& t, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    TPerTemplateResultsMap::iterator it = m_PerTemplateResults.find( t.GetName() );
    if( it == m_PerTemplateResults.end() )
    {
        m_PerTemplateResults[ t.GetName() ] = TResultsMap();
    }

    StoreResult( m_PerTemplateResults[ t.GetName() ], path, name, result );
}

void COperatorResultStore::StoreGlobalResult( const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    StoreResult( m_GlobalResults, path, name, result );
}

std::weak_ptr<const COperatorResult> COperatorResultStore::GetPerTemplateResult( const CTemplate& t, const std::string& name ) const
{
    TPerTemplateResultsMap::const_iterator it1 = m_PerTemplateResults.find( t.GetName() );
    if( it1 != m_PerTemplateResults.end() )
    {
        const TResultsMap& m = (*it1).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            return (*it2).second;
        }
    }
    
    Log.LogError("COperatorResultStore::GetPerTemplateResult, per template result %s not found",name.c_str());
    //throw runtime_error("COperatorResultStore::GetPerTemplateResult, global result not found");
    return std::weak_ptr<const COperatorResult>();
}

TOperatorResultMap COperatorResultStore::GetPerTemplateResult( const std::string& name ) const
{
    TOperatorResultMap map;
    TPerTemplateResultsMap::const_iterator it;
    for( it = m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); ++it )
    {
        std::string tplName = (*it).first;

        const TResultsMap& m = (*it).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            map[tplName] = (*it2).second;
        }
    }

    return map;
}

/**
 * /brief Returns the best global result, if there is one.
 * 
 * The somewhat strange syntax of this method is due to the usage of std::map, which will yield an iterator when accessing a member using .find().
 * This method will find the global result entry with key "name", and if it exists, will return its .second member (which will be a COperatorResult or derived object). Otherwise, will return a pointer to an empty COperatorResult.
 */

std::weak_ptr<const COperatorResult> COperatorResultStore::GetGlobalResult( const std::string& name ) const
{
    TResultsMap::const_iterator it = m_GlobalResults.find( name );
    if( it != m_GlobalResults.end() )
    {
      return (*it).second;
    }
    else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter()<<"Unknown global result:"<<name);
    //Log.LogError("COperatorResultStore::GetGlobalResult, global result %s not found",name.c_str());
    //throw runtime_error("COperatorResultStore::GetGlobalResult, global result not found");

}

std::weak_ptr<const COperatorResult> COperatorResultStore::GetScopedGlobalResult( const std::string& name ) const
{
  return GetGlobalResult(GetScopedName(name));
}

std::weak_ptr<const COperatorResult> COperatorResultStore::GetGlobalResult( const std::string& objectType,
                                                                            const std::string& method,
                                                                            const std::string& name ) const
{
    std::ostringstream oss;
    oss << objectType << "." << method << "." << name;
    return GetGlobalResult(oss.str());
}


std::weak_ptr<const COperatorResult> COperatorResultStore::GetGlobalResult( const std::string& objectType,
                                                                            const std::string& method,
                                                                            const std::string& name ,
                                                                            const int& rank) const
{
    std::ostringstream oss;
    oss << objectType << "." << method << "." << name<<rank;
    return GetGlobalResult(oss.str());
}



/**
 * @brief COperatorResultStore::CreateResultStorage
 *
 * @return 0 if already exists, 1 if just created, -1 if error
 */
Int32 COperatorResultStore::CreateResultStorage( std::fstream& stream, const bfs::path& path, const bfs::path& baseDir ) const
{
    Int32 ret = -1;

    bfs::path outputFilePath = bfs::path( baseDir );
    outputFilePath /= path.string();


    if( bfs::exists( outputFilePath.parent_path() ) == false )
    {
        bfs::create_directories( outputFilePath.parent_path() );
    }

    if( bfs::exists( outputFilePath ) == false )
    {
        ret = 1;
    }else{
        ret = 0;
    }

    stream.open( outputFilePath.string().c_str(), std::fstream::out | std::fstream::app);
    if( stream.rdstate() & std::ios_base::failbit )
    {
        return -1;
    }

    return ret;
}

void COperatorResultStore::SaveRedshiftResult( const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "redshift.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tsnrOII\tlfOII\tType"<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "redshiftresult" ).lock();
        if(result){
            result->SaveLine( outputStream );
        }else{
            throw std::runtime_error("Unable to retrieve redshift result for saving");
        }
    }
}

void COperatorResultStore::SaveRedshiftResultError(  const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "redshift.csv" ), dir );

        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tsnrOII\tlfOII\tType"<< std::endl;
        }


        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1"<< std::endl;
    }
}

void COperatorResultStore::SaveStellarResult(const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "stellar.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tsnrOII\tlfOII\tType"<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "stellarsolve.stellarresult" ).lock();
        if(result){
            result->SaveLine(outputStream );
        }//else{
        //    throw std::runtime_error("Unable to retrieve stellar result for saving");
        //}
    }
}

void COperatorResultStore::SaveStellarResultError(  const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "stellar.csv" ), dir );

        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tsnrOII\tlfOII\tType"<< std::endl;
        }


        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1"<< std::endl;
    }
}


void COperatorResultStore::SaveQsoResult( const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "qso.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tsnrOII\tlfOII\tType"<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "qsosolve.qsoresult" ).lock();
        if(result){
            result->SaveLine(outputStream );
        }//else{
        //    throw std::runtime_error("Unable to retrieve stellar result for saving");
        //}
    }
}

void COperatorResultStore::SaveQsoResultError(  const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "qso.csv" ), dir );

        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift\tMerit\tTemplate\tMethod\tDeltaz\tReliability\tsnrHa\tlfHa\tsnrOII\tlfOII\tType"<< std::endl;
        }

        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t-1"<< std::endl;
    }
}


void COperatorResultStore::SaveClassificationResult(const bfs::path& dir )
{
    // Append classif. result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "classification.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\ttype\tLogEvidenceG\tLogEvidenceS\tLogEvidenceQ\tProbaG\tProbaS\tProbaQ"<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "classificationresult" ).lock();
        if(result){
            result->SaveLine( outputStream );
        }//else{
        //    throw std::runtime_error("Unable to retrieve classification result for saving");
        //}
    }
}

void COperatorResultStore::SaveClassificationResultError(  const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append best classification result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "classification.csv" ), dir );

        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\ttype\tevidenceG\tevidenceS\tevidenceQ"<< std::endl;
        }


        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1"<< std::endl;
    }
}

void COperatorResultStore::SaveCandidatesResult(  const bfs::path& dir )
{
    // Append candidate result line to output candidate file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "candidates.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRank_1\tID_1\tRedshift_1\tProb_1\tRank_PDF_1\tgaussAmp_1\tgaussSigma_1\tRank_2\tID_2\tRedshift_2\tProb_2\tRank_PDF_2\tgaussAmp_2\tgaussSigma_2\t..."<< std::endl;
        }
        //*/

        auto  result = GetGlobalResult( "candidatesresult" ).lock();
        if(result){
            result->SaveLine( outputStream );
        }
    }
}

void COperatorResultStore::SaveCandidatesResultError( const std::string spcName,  const std::string processingID, const bfs::path& dir )
{
    // Append candidate result line to output candidate file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        Int32 ret = CreateResultStorage( outputStream, bfs::path( "candidates.csv" ), dir );

        //*
        if(ret==1)
        {
            outputStream <<  "#Spectrum\tProcessingID\tRedshift_1\tProb_1\tgaussAmp_1\tgaussSigma_1\tRedshift_2\tProb_2\tgaussAmp_2\tgaussSigma_2\t..."<< std::endl;
        }
        //*/

        outputStream <<  spcName << "\t" << processingID << "\t-1\t-1\t-1\t-1\t-1\t-1\t-1"<< std::endl;
    }
}

void COperatorResultStore::SaveReliabilityResult( const bfs::path& dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        CreateResultStorage( outputStream, bfs::path( "redshiftreliability.csv" ), dir );

        auto  result = GetGlobalResult( "zReliability/result.zpredict" ).lock();
        if(result){
            result->SaveLine( outputStream );
        }
    }
}


void COperatorResultStore::SaveAllResults(const bfs::path& dir, const std::string opt ) const
{
    std::string opt_lower = opt;
    boost::algorithm::to_lower(opt_lower);

    // Store global result
    {
        TResultsMap::const_iterator it;
        for( it=m_GlobalResults.begin(); it != m_GlobalResults.end(); it++ )
        {
            std::string resultName = (*it).first;
            Log.LogInfo("save result %s",resultName.c_str()); 
            auto  result = (*it).second;

            bool firstpass = false;
            bool saveThisResult = false;
	        bool saveJSON = false;
            
            if(opt_lower=="all" || opt_lower=="global"){
                saveThisResult = true;
                //save extrema results
                std::string extremaresTagRes = "linemodel_extrema";
                std::size_t foundstr = resultName.find(extremaresTagRes.c_str());
                if (foundstr!=std::string::npos){
                  Log.LogInfo("save only json");                   
                  saveThisResult = false;
                  saveJSON = true;
                }
            }
            else if(opt_lower=="linemeas")
            {
                std::string linemeasTagRes = "linemodel_fit_extrema_0";
                std::size_t foundstr = resultName.find(linemeasTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
            }else if(opt_lower=="default" || opt_lower=="all" || opt_lower=="global")
            {
                //save extrema results
                std::string extremaresTagRes = "linemodel_extrema";
                std::size_t foundstr = resultName.find(extremaresTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=false;
		            saveJSON = true;
                }

                //save first pass extrema results
                std::string firstpassextremaresTagRes = "linemodel_firstpass_extrema";
                foundstr = resultName.find(firstpassextremaresTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }

                //save best fit parameters
                std::string linefitTagRes = "linemodel_fit_extrema_0";
                foundstr = resultName.find(linefitTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
                //save best model
                std::string linemodelTagRes = "linemodel_spc_extrema_0";
                foundstr = resultName.find(linemodelTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
                //save best continuum
                std::string continuummodelTagRes = "linemodel_continuum_extrema_0";
                foundstr = resultName.find(continuummodelTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
                //save redshiftresult brief
                std::string redshiftresultTagRes = "redshiftresult";
                foundstr = resultName.find(redshiftresultTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }

                //save pdf
                std::string pdfTagRes = "logposterior.logMargP_Z_data";
                foundstr = resultName.find(pdfTagRes.c_str());
                if (foundstr!=std::string::npos){
                    saveThisResult=true;
                }
            }

            if (resultName == "linemodel_extrema" ||  // do not save linemodel_extrema.csv, json works well
                resultName == "linemodelsolve.linemodel" || // contains useleess chi2 and duplicate info from linemodel_extrema 
                resultName.find("linemodel_rules_extrema") !=std::string::npos ) // empty files...
              {
                Log.LogInfo("do not save csv");
                saveThisResult=false;
              
              }
              /*            
            if(!saveThisResult && !saveJSON)
            {
                continue;
            }


            std::string firstpassextremaresTagRes = "linemodel_firstpass_extrema";
            std::size_t foundstr = resultName.find(firstpassextremaresTagRes.c_str());
            if (foundstr!=std::string::npos){
                    firstpass=true;
            }
              */
            std::fstream outputStream;
            // Save result at root of output directory

            if (saveThisResult)
              {
                CreateResultStorage( outputStream, bfs::path( resultName + ".csv"), dir );
                result->Save(outputStream);
              }
	        if(saveJSON){
		        std::fstream outputJSONStream;
                TResultsMap::const_iterator it = m_GlobalResults.find("candidatesresult" );
 
		        CreateResultStorage( outputJSONStream, bfs::path( resultName + ".json"), dir );
		        result->SaveJSON(outputJSONStream);
	        }
        }
    }

    // Store per template results
    // Do not store per templates result
    /*
    if(opt_lower=="all"){
        TPerTemplateResultsMap::const_iterator it;
        for( it=m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); it++ )
        {
            std::string templateName = (*it).first;
            const TResultsMap& resultMap = (*it).second;

            TResultsMap::const_iterator it2;
            for( it2=resultMap.begin(); it2 != resultMap.end(); it2++ )
            {
                std::string resultName = (*it2).first;
                auto  result = (*it2).second;

                std::fstream outputStream;
                // Save result in sub directories of output directory
                bfs::path outputFilePath = bfs::path( templateName );
                outputFilePath /= std::string( resultName + ".csv" );
                CreateResultStorage( outputStream, outputFilePath, dir );
                result->Save( store, outputStream );
            }
        }
    }
    */
}


std::string COperatorResultStore::GetScope( const COperatorResult&  result) const
{
    std::string n="";

    TResultsMap::const_iterator it;
    for( it=m_GlobalResults.begin(); it != m_GlobalResults.end(); it++ )
    {
        auto  r = (*it).second;
        if(&result == r.get() ){
            std::string s = (*it).first;

            std::size_t found = s.rfind(".");
            if (found!=std::string::npos){
                n = s.substr(0,found);
                n = n.append(".");
            }
            break;
        }
    }

    return n;
}


void COperatorResultStore::getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, Float64& v) const
{
 std:weak_ptr<const COperatorResult> result;
    if (method == "linemodelsolve")
    {
      result=GetGlobalResult(object_type,method,"linemodel_extrema");
      return result.lock()->getCandidateData(rank,name,v);
    }
  else if(method.compare("templatefittingsolve") == 0)
    {
      if (name == "Redshift" || name == "RedshiftError" || name =="RedshiftProba")
        {
          result = GetGlobalResult(object_type,method,"candidatesresult");
          return result.lock()->getCandidateData(rank,name,v);
        }
      else result=GetGlobalResult(object_type,method,"templatefitting_fitcontinuum_extrema_",rank);      
    }
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown object_type "<<object_type);
  result.lock()->getData(name,v);
}

void COperatorResultStore::getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, std::string& v) const
{
 std:weak_ptr<const COperatorResult> result;
  if (method == "linemodelsolve")
    {
      result=GetGlobalResult(object_type,method,"linemodel_extrema");
      return result.lock()->getCandidateData(rank,name,v);
    }
  else if(method.compare("templatefittingsolve") == 0) result=GetGlobalResult(object_type,method,"templatefitting_fitcontinuum_extrema_",rank);
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown object_type "<<object_type);
  result.lock()->getData(name,v);
}

void COperatorResultStore::getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, Int32& v) const
{
 std:weak_ptr<const COperatorResult> result;
  if (name == "Rank") result= GetGlobalResult(object_type,method,"candidatesresult");
  else if (method == "linemodelsolve")
    {
      result=GetGlobalResult(object_type,method,"linemodel_extrema");
      return result.lock()->getCandidateData(rank,name,v);
    }
  else if(method.compare("templatefittingsolve") == 0) result=GetGlobalResult(object_type,method,"templatefitting_fitcontinuum_extrema_",rank);
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown object_type "<<object_type);
  if (name == "Rank") result.lock()->getCandidateData(rank,name,v);
  else result.lock()->getData(name,v);
}

void COperatorResultStore::getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, double **data, int *size) const
{
 std:weak_ptr<const COperatorResult> result;
  if (name.find("Model") != std::string::npos)
    {
      if (method == "templatefittingsolve") result = GetGlobalResult(object_type,method,"templatefitting_spc_extrema_",rank);
      else if(method == "linemodelsolve") result = GetGlobalResult(object_type,method,"linemodel_spc_extrema_",rank);
      else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown method "<<method<< " for attribute" <<name );

    }
  else if (name.find("FittedRays") != std::string::npos) result = GetGlobalResult(object_type,method,  "linemodel_fit_extrema_", rank); 
  else if (name.find("BestContinuum") != std::string::npos)  result = GetGlobalResult(object_type,method,"linemodel_continuum_extrema_", rank);
  else if (name.compare("ContinuumIndexesColor") == 0 || name.compare("ContinuumIndexesBreak") == 0)
    {
      result = GetGlobalResult(object_type,method,"linemodel_extrema");
      result.lock()->getCandidateData(rank,name,data,size);
      return;
    }
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown data "<<name);
    
  result.lock()->getData(name,data,size);
}

/*
void COperatorResultStore::getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, std::string *data, int *size) const
{
 std:weak_ptr<const COperatorResult> result;
  std::ostringstream oss;
  if (name.find("model_") != std::string::npos)  oss << "linemodelsolve.linemodel_spc_extrema_"<< rank;
  else if (name.find("fitted_rays_") != std::string::npos)  oss << "linemodelsolve.linemodel_fit_extrema_"<< rank;
  else
    {
      Log.LogError("unknown data "<<name);
    }
  result = GetGlobalResult(oss.str());
  result.lock()->getData(name,data,size);
}
*/


void COperatorResultStore::getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, int **data, int *size) const 
{
 std:weak_ptr<const COperatorResult> result; 

  if (name.find("Model") != std::string::npos) result = GetGlobalResult(object_type,method,"templatefitting_spc_extrema_",rank);
  else if (name.find("FittedRays") != std::string::npos) result = GetGlobalResult(object_type,method,  "linemodel_fit_extrema_", rank);
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown data "<<name);
    
  result.lock()->getData(name,data,size);
}


void COperatorResultStore::getData(const std::string& object_type,const std::string& method,const std::string& name, Int32& v) const
{
 std:weak_ptr<const COperatorResult> result = GetGlobalResult(object_type,method,"candidatesresult");
  if (name == "NbCandidates") result.lock()->getData(name,v);
  else throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() <<"unknown object_type "<<object_type);
      
}

void COperatorResultStore::getData(const std::string& object_type,const std::string& method,const std::string& name, Float64& v) const
{
  std:weak_ptr<const COperatorResult> result;
  if(name.compare("snrHa") == 0 || name.compare("lfHa") == 0 ||
     name.compare("snrOII") == 0 || name.compare("lfOII") == 0) result = GetGlobalResult("redshiftresult");
  else if(object_type=="classification") result = GetGlobalResult("classificationresult");
  else result = GetGlobalResult(object_type,method,"candidatesresult");
  result.lock()->getData(name,v);
}

void COperatorResultStore::getData(const std::string& object_type,const std::string& method,const std::string& name, std::string& v) const
{
  std:weak_ptr<const COperatorResult> result;
  if (object_type.compare("classification") == 0) result = GetGlobalResult("classificationresult");
  else result = GetGlobalResult(object_type,method,"candidatesresult");
  result.lock()->getData(name,v);
}

void COperatorResultStore::getData(const std::string& object_type,const std::string& method,const std::string& name,double **data, int *size) const
{
  std:weak_ptr<const COperatorResult> result = GetGlobalResult(object_type,method,"pdf");
  result.lock()->getData(name,data,size);
}

void COperatorResultStore::test()
{
  std::shared_ptr<CPdfMargZLogResult> testResult = std::shared_ptr<CPdfMargZLogResult>(new CPdfMargZLogResult());
  testResult->Redshifts.clear();
  testResult->Redshifts.push_back(1.2);
  testResult->Redshifts.push_back(3.4);
  testResult->Redshifts.push_back(7.9);

  testResult->valProbaLog.clear();
  testResult->valProbaLog.push_back(27.4);
  testResult->valProbaLog.push_back(67.8);
  testResult->valProbaLog.push_back(92.7);
  
  StoreGlobalResult("","zPDF/logposterior.logMargP_Z_data",testResult);

  
}


void  COperatorResultStore::StoreScopedPerTemplateResult( const CTemplate& t, const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    StorePerTemplateResult( t, GetCurrentScopeName(), name, result );
}

void COperatorResultStore::StoreScopedGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    StoreGlobalResult( GetCurrentScopeName(), name, result );
}

void COperatorResultStore::StoreGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult> result )
{
    StoreGlobalResult( "", name, result );
}



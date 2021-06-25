#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/statistics/pdfcandidateszresult.h"
#include "RedshiftLibrary/method/classificationresult.h"
#include "RedshiftLibrary/linemodel/linemodelextremaresult.h"
#include "RedshiftLibrary/operator/extremaresult.h"
#include "RedshiftLibrary/operator/tplCombinationExtremaResult.h"
#include "RedshiftLibrary/linemodel/modelfittingresult.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/spectraFluxResult.h"


#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"

using namespace NSEpic;

namespace bfs = boost::filesystem;


COperatorResultStore::COperatorResultStore(const TScopeStack& scope):
  CScopeStore(scope)
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


std::shared_ptr<const CClassificationResult> COperatorResultStore::GetClassificationResult( const std::string& objectType,
                                                                            const std::string& method,
                                                                            const std::string& name ) const
{
    std::ostringstream oss;
    oss << objectType << "." << method << "." << name;
    std::weak_ptr<const COperatorResult> cor = GetGlobalResult(oss.str());

    return std::dynamic_pointer_cast<const CClassificationResult>(GetGlobalResult(oss.str()).lock());
}

std::shared_ptr<const CPdfMargZLogResult> COperatorResultStore::GetPdfMargZLogResult( const std::string& objectType,
                                                                            const std::string& method,
                                                                            const std::string& name ) const
{
    std::ostringstream oss;
    oss << objectType << "." << method << "." << name;
    std::weak_ptr<const COperatorResult> cor = GetGlobalResult(oss.str());

    return std::dynamic_pointer_cast<const CPdfMargZLogResult>(GetGlobalResult(oss.str()).lock());
}

//std::shared_ptr<const CLineModelExtremaResult<TLineModelResult>>

//TODO AA those getters should have an argument std::string dataset, rather than hard coded values

std::shared_ptr<const TLineModelResult> COperatorResultStore::GetLineModelResult(const std::string& objectType,
										 const std::string& method,
										 const std::string& name ,
										 const int& rank
										 ) const
    
{
  std::shared_ptr<const COperatorResult> cop = GetGlobalResult(objectType,method,name).lock()->getCandidate(rank,"model_parameters");
  std::shared_ptr<const TLineModelResult> tlm = std::dynamic_pointer_cast<const TLineModelResult>(cop); 
    return tlm;



}

std::shared_ptr<const TTplCombinationResult> COperatorResultStore::GetTplCombinationResult(const std::string& objectType,
										 const std::string& method,
										 const std::string& name ,
										 const int& rank
										 ) const
    
{
  std::shared_ptr<const COperatorResult> cop = GetGlobalResult(objectType,method,name).lock()->getCandidate(rank,"model_parameters");
  std::shared_ptr<const TTplCombinationResult> ttc = std::dynamic_pointer_cast<const TTplCombinationResult>(cop); 
  return ttc;
}

std::shared_ptr<const TExtremaResult> COperatorResultStore::GetExtremaResult(const std::string& objectType,
										 const std::string& method,
										 const std::string& name ,
										 const int& rank
										 ) const
    
{
  std::shared_ptr<const COperatorResult> cop = GetGlobalResult(objectType,method,name).lock()->getCandidate(rank,"model_parameters");
  std::shared_ptr<const TExtremaResult> tlm = std::dynamic_pointer_cast<const TExtremaResult>(cop); 
    return tlm;
}

std::shared_ptr<const CModelFittingResult> COperatorResultStore::GetModelFittingResult(const std::string& objectType,
										       const std::string& method,
										       const std::string& name ,
										       const int& rank
										       ) const
    
{
  return std::dynamic_pointer_cast<const CModelFittingResult>(GetGlobalResult(objectType,method,name).lock()->getCandidate(rank,"fitted_rays"));
}

std::shared_ptr<const CModelFittingResult> COperatorResultStore::GetModelFittingResult(const std::string& objectType,
										       const std::string& method,
										       const std::string& name 
										       ) const
    
{
  return std::dynamic_pointer_cast<const CModelFittingResult>(GetGlobalResult(objectType,method,name).lock());
}

std::shared_ptr<const CSpectraFluxResult> COperatorResultStore::GetSpectraFluxResult(const std::string& objectType,
										     const std::string& method,
										     const std::string& name ,
										     const int& rank
										     ) const
{
  return std::dynamic_pointer_cast<const CSpectraFluxResult>(GetGlobalResult(objectType,method,name).lock()->getCandidate(rank,"continuum"));
}

std::shared_ptr<const CModelSpectrumResult> COperatorResultStore::GetModelSpectrumResult(const std::string& objectType,
										       const std::string& method,
										       const std::string& name ,
										       const int& rank
										       ) const
    
{
  return std::dynamic_pointer_cast<const CModelSpectrumResult>(GetGlobalResult(objectType,method,name).lock()->getCandidate(rank,"model"));
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

const std::string&  COperatorResultStore::GetGlobalResultType(const std::string& objectType,
                                                              const std::string& method,
                                                              const std::string& name ) const
{
  return GetGlobalResult(objectType,method,name).lock()->getType();
}

const std::string&  COperatorResultStore::GetCandidateResultType(const std::string& objectType,
                                                              const std::string& method,
								 const std::string& name ,
								 const std::string& dataset) const
{
  return GetGlobalResult(objectType,method,name).lock()->getCandidateDatasetType(dataset);
}

bool COperatorResultStore::HasCandidateDataset(const std::string& objectType,
						       const std::string& method,
						       const std::string& name ,
						       const std::string& dataset) const
{
  if (HasDataset(objectType,method,name))
    {
      return  GetGlobalResult(objectType,method,name).lock()->HasCandidateDataset(dataset);
    }
  else return false;
}

bool COperatorResultStore::HasDataset(const std::string& objectType,
                                      const std::string& method,
                                      const std::string& name ) const
{
  std::ostringstream oss;
  oss << objectType << "." << method << "." << name;
  TResultsMap::const_iterator it = m_GlobalResults.find( oss.str() );
  return (it != m_GlobalResults.end());
}



int COperatorResultStore::getNbRedshiftCandidates(const std::string& objectType,
						  const std::string& method) const
{
  std::ostringstream oss;
  oss << objectType << "." << method << ".candidatesresult" ;
  TResultsMap::const_iterator it = m_GlobalResults.find( oss.str() );
  if (it != m_GlobalResults.end())
    {
      std::weak_ptr<const COperatorResult> cor = GetGlobalResult(oss.str());

      return std::dynamic_pointer_cast<const PdfCandidatesZResult>(GetGlobalResult(oss.str()).lock())->size();
    }
  else return 0;

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
  
  StoreGlobalResult("","star.templatefittingsolve.solveResult",testResult);

  std::shared_ptr<const COperatorResult> cor = GetSolveResult("star");
  
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


std::vector<std::string> COperatorResultStore::getProcessedObjectTypes() const
{
  std::vector<std::string> res;

  return res;
}

std::shared_ptr<const COperatorResult> COperatorResultStore::GetSolveResult(const std::string& objectType) const
{
      TResultsMap::const_iterator it;
    for( it=m_GlobalResults.begin(); it != m_GlobalResults.end(); it++ )
    {
      std::string s = (*it).first;
      std::size_t first_dot = s.find(".");
      std::size_t second_dot = s.rfind(".");
      std::string object_type = s.substr(0,first_dot);
      std::string rs_name = s.substr(second_dot+1,s.size()-second_dot-1);
      if(object_type == objectType && rs_name == "solveResult") 
        return it->second;
    }
    throw GlobalException(UNKNOWN_ATTRIBUTE,Formatter() << "no solveResult found for objectType" <<objectType);
}

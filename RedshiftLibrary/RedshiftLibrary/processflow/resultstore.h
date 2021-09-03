// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_PROCESSFLOW_OPERATORRESULTSTORE_
#define _REDSHIFT_PROCESSFLOW_OPERATORRESULTSTORE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/processflow/scopestore.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <vector>
#include <ostream>

namespace NSEpic
{

class CTemplate;
class CClassificationResult;
class CReliabilityResult;
class CPdfMargZLogResult;
class TLineModelResult;
class TTplCombinationResult;
class TExtremaResult;
class CModelFittingResult;
class CModelSpectrumResult;
class CSpectraFluxResult;
class CLineModelSolution;
template<class T=TLineModelResult>  class CLineModelExtremaResult;
template<class T=TTplCombinationResult>  class CTplCombinationExtremaResult;
  //  class CLineModelExtremaResult<TLineModelResult>;
  //class LineModelExtremaResult;
  /**
 * \ingroup Redshift
 */
class COperatorResultStore : public CScopeStore
{

public:

    typedef std::map< std::string, std::shared_ptr< const COperatorResult> >    TResultsMap;
    typedef std::map< std::string, TResultsMap>                                 TPerTemplateResultsMap;


    COperatorResultStore(const TScopeStack& scopeStack);

    void                    StorePerTemplateResult( const CTemplate& t, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );
    void                    StoreGlobalResult( const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );

    std::weak_ptr<const COperatorResult> GetPerTemplateResult( const CTemplate& t, const std::string& name ) const;
    TOperatorResultMap      GetPerTemplateResult( const std::string& name ) const;
    std::weak_ptr<const COperatorResult>  GetGlobalResult( const std::string& name ) const;
    std::weak_ptr<const COperatorResult>  GetScopedGlobalResult( const std::string& name ) const;

  std::weak_ptr<const COperatorResult>  GetGlobalResult(const std::string& objectType,
                                                        const std::string& method,
                                                        const std::string& name ) const;
  std::weak_ptr<const COperatorResult>  GetGlobalResult(const std::string& objectType,
                                                        const std::string& method,
                                                        const std::string& name ,
                                                        const int& rank) const;

  const std::string&  GetGlobalResultType(const std::string& objectType,
                                          const std::string& method,
                                          const std::string& name ) const;

  const std::string&  GetCandidateResultType(const std::string& objectType,
					     const std::string& method,
					     const std::string& name ,
					     const std::string& dataset) const;

  bool HasCandidateDataset(const std::string& objectType,
			   const std::string& method,
			   const std::string& name ,
			   const std::string& dataset) const;

  bool HasDataset(const std::string& objectType,
                  const std::string& method,
                  const std::string& name ) const;

  std::shared_ptr<const CClassificationResult> GetClassificationResult(const std::string& objectType,
                                                                            const std::string& method,
                                                                            const std::string& name ) const;
    std::shared_ptr<const CReliabilityResult> GetReliabilityResult(const std::string& objectType,
                                                                            const std::string& method,
                                                                            const std::string& name ) const;

  std::shared_ptr<const CPdfMargZLogResult> GetPdfMargZLogResult(const std::string& objectType,
								    const std::string& method,
								    const std::string& name ) const;

  std::shared_ptr<const TLineModelResult> GetLineModelResult(const std::string& objectType,
							     const std::string& method,
							     const std::string& name ,
							     const int& rank
							     ) const;
  std::shared_ptr<const TTplCombinationResult> GetTplCombinationResult(const std::string& objectType,
										 const std::string& method,
										 const std::string& name ,
										 const int& rank
										 ) const;
  std::shared_ptr<const TExtremaResult> GetExtremaResult(const std::string& objectType,
										 const std::string& method,
										 const std::string& name ,
										 const int& rank
									       ) const;

  std::shared_ptr<const CModelFittingResult> GetModelFittingResult(const std::string& objectType,
								   const std::string& method,
								   const std::string& name ,
								   const int& rank
								   ) const  ;

  std::shared_ptr<const CModelFittingResult> GetModelFittingResult(const std::string& objectType,
								   const std::string& method,
								   const std::string& name 
								   ) const  ;
  
  std::shared_ptr<const CLineModelSolution> GetLineModelSolution(const std::string& objectType,
								   const std::string& method,
								   const std::string& name 
								   ) const  ;

  std::shared_ptr<const CModelSpectrumResult> GetModelSpectrumResult(const std::string& objectType,
								     const std::string& method,
								     const std::string& name ,
								     const int& rank
								     ) const  ;
  std::shared_ptr<const CModelSpectrumResult> GetModelSpectrumResult(const std::string& objectType,
								     const std::string& method,
								     const std::string& name 
								     ) const  ;

  std::shared_ptr<const CSpectraFluxResult> GetSpectraFluxResult(const std::string& objectType,
								   const std::string& method,
								   const std::string& name ,
								   const int& rank
								   ) const  ;

  std::vector<std::string> getProcessedObjectTypes() const;
  
  int getNbRedshiftCandidates(const std::string& objectType,const std::string& method) const;

  std::shared_ptr<const COperatorResult> GetSolveResult(const std::string& objectType) const;
  
  
  //From DataStore, above should be removed and integrated into these
      void                    StoreGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult> result );

    void                            StoreScopedPerTemplateResult( const CTemplate& t, const std::string& name, std::shared_ptr<const COperatorResult>  result );
    void                            StoreScopedGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult>  result );
  
  
    std::string             GetScope( const COperatorResult&  result) const;
  
  void test();

  
protected:

    void StoreResult( TResultsMap& map, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );

    Int32 CreateResultStorage( std::fstream& stream, const boost::filesystem::path& path, const boost::filesystem::path& baseDir ) const;

    TPerTemplateResultsMap          m_PerTemplateResults;
    TResultsMap                     m_GlobalResults;

};


}

#endif

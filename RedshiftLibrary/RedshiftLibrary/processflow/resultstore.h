#ifndef _REDSHIFT_OPERATOR_RESULTSTORE_
#define _REDSHIFT_OPERATOR_RESULTSTORE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/exception.h>
#include <RedshiftLibrary/processflow/result.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <vector>
#include <ostream>

namespace NSEpic
{

class CTemplate;
class CDataStore;

/**
 * \ingroup Redshift
 */
class COperatorResultStore
{

public:

    typedef std::map< std::string, std::shared_ptr< const COperatorResult> >    TResultsMap;
    typedef std::map< std::string, TResultsMap>                                 TPerTemplateResultsMap;


    COperatorResultStore();
    virtual ~COperatorResultStore();

    void                    StorePerTemplateResult( const CTemplate& t, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );
    void                    StoreGlobalResult( const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );

    std::weak_ptr<const COperatorResult> GetPerTemplateResult( const CTemplate& t, const std::string& name ) const;
    TOperatorResultMap      GetPerTemplateResult( const std::string& name ) const;
    std::weak_ptr<const COperatorResult>  GetGlobalResult( const std::string& name ) const;
    void                    DeleteGlobalResult(const std::string& path, const std::string& name);

    void                    SaveRedshiftResultError( const std::string spcName, const std::string processingID, const boost::filesystem::path& dir );
    void                    SaveRedshiftResult( const CDataStore& store, const boost::filesystem::path& dir );
    void                    SaveCandidatesResult( const CDataStore& store, const boost::filesystem::path& dir );
    void                    SaveCandidatesResultError( const std::string spcName, const std::string processingID, const boost::filesystem::path& dir );
    void                    SaveAllResults(const CDataStore& store, const boost::filesystem::path& dir , const std::string opt) const;
    void                    SaveReliabilityResult( const CDataStore& store, const boost::filesystem::path& dir );
    void                    SaveStellarResultError( const std::string spcName, const std::string processingID, const boost::filesystem::path& dir );
    void                    SaveStellarResult( const CDataStore& store, const boost::filesystem::path& dir );
    void                    SaveQsoResultError( const std::string spcName, const std::string processingID, const boost::filesystem::path& dir );
    void                    SaveQsoResult( const CDataStore& store, const boost::filesystem::path& dir );
    void                    SaveClassificationResultError( const std::string spcName, const std::string processingID, const boost::filesystem::path& dir );
    void                    SaveClassificationResult( const CDataStore& store, const boost::filesystem::path& dir );

    std::string             GetScope( const COperatorResult&  result) const;


  // it should be get<Type><Method>CandidateParam(int rank, std::string param, std::string method)
  void getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, Float64& v) const;
  void getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, Int32& v) const;
  void getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, std::string& v) const;
  void getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, double **data, int *size) const;
  //void getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, std::string *data, int *size) const;
  void getCandidateData(const std::string& object_type,const std::string& method,const int& rank,const std::string& name, int **data, int *size) const;

  void getData(const std::string& object_type,const std::string& method,const std::string& name, Int32& v) const;
  void getData(const std::string& object_type,const std::string& method,const std::string& name, Float64& v) const;
  void getData(const std::string& object_type,const std::string& method,const std::string& name, std::string& v) const;
  void getData(const std::string& object_type,const std::string& method,const std::string& name, double **data, int *size) const;
  void getData(const std::string& object_type,const std::string& method,const std::string& name, int **data, int *size) const;

  void test();
protected:

    void StoreResult( TResultsMap& map, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );

    Int32 CreateResultStorage( std::fstream& stream, const boost::filesystem::path& path, const boost::filesystem::path& baseDir ) const;

    TPerTemplateResultsMap          m_PerTemplateResults;
    TResultsMap                     m_GlobalResults;

};


}

#endif

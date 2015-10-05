#ifndef _REDSHIFT_OPERATOR_RESULTSTORE_
#define _REDSHIFT_OPERATOR_RESULTSTORE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/processflow/result.h>

#include <boost/filesystem.hpp>
#include <vector>
#include <ostream>

namespace NSEpic
{

class CTemplate;
class CDataStore;

class COperatorResultStore
{

public:

    typedef std::map< std::string, CConstRef<COperatorResult> > TResultsMap;
    typedef std::map< std::string, TResultsMap>                 TPerTemplateResultsMap;


    COperatorResultStore();
    virtual ~COperatorResultStore();

    Void                    StorePerTemplateResult( const CTemplate& t, const char* path, const char* name, const COperatorResult& result );
    Void                    StoreGlobalResult( const char* path, const char* name, const COperatorResult& result );

    const COperatorResult*  GetPerTemplateResult( const CTemplate& t, const char* name ) const;
    TOperatorResultMap      GetPerTemplateResult( const char* name ) const;
    const COperatorResult*  GetGlobalResult( const char* name ) const;

    Void                    SaveRedshiftResultHeader( const char* dir );
    Void                    SaveRedshiftResult( const CDataStore& store, const char* dir );
    Void                    SaveAllResults( const CDataStore& store, const char* dir ) const;

    std::string             GetScope(CConstRef<COperatorResult>  result) const;

protected:

    void StoreResult( TResultsMap& map, const char* path, const char* name, const COperatorResult& result );

    void CreateResultStorage( std::fstream& stream, const boost::filesystem::path& path, const boost::filesystem::path& baseDir ) const;

    TPerTemplateResultsMap          m_PerTemplateResults;
    TResultsMap                     m_GlobalResults;
    std::string                     m_SpectrumName;

};


}

#endif

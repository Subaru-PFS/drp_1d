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

    Void                    StorePerTemplateResult( const CTemplate& t, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );
    Void                    StoreGlobalResult( const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );

    std::weak_ptr<const COperatorResult> GetPerTemplateResult( const CTemplate& t, const std::string& name ) const;
    TOperatorResultMap      GetPerTemplateResult( const std::string& name ) const;
    std::weak_ptr<const COperatorResult>  GetGlobalResult( const std::string& name ) const;

    Void                    SaveRedshiftResultHeader( const boost::filesystem::path& dir );
    Void                    SaveRedshiftResult( const CDataStore& store, const boost::filesystem::path& dir );
    Void                    SaveAllResults( const CDataStore& store, const boost::filesystem::path& dir ) const;

    std::string             GetScope( const COperatorResult&  result) const;

protected:

    void StoreResult( TResultsMap& map, const std::string& path, const std::string& name, std::shared_ptr<const COperatorResult> result );

    void CreateResultStorage( std::fstream& stream, const boost::filesystem::path& path, const boost::filesystem::path& baseDir ) const;

    TPerTemplateResultsMap          m_PerTemplateResults;
    TResultsMap                     m_GlobalResults;

};


}

#endif

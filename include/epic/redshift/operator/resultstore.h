#ifndef _REDSHIFT_OPERATOR_RESULTSTORE_
#define _REDSHIFT_OPERATOR_RESULTSTORE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>
#include <epic/redshift/operator/result.h>

#include <boost/filesystem.hpp>
#include <vector>
#include <ostream>

namespace NSEpic
{

class CTemplate;

class COperatorResultStore : public CManagedObject
{

public:

    typedef std::map< std::string, CConstRef<COperatorResult> > TResultsMap;
    typedef std::map< std::string, TResultsMap>                 TPerTemplateResultsMap;
    typedef std::vector<std::string>                            TScopeStack;

    class CAutoScope {
    public:
        CAutoScope( COperatorResultStore& store, const char* name );
        ~CAutoScope();
    private:
        COperatorResultStore* m_Store;
    };

    COperatorResultStore();
    virtual ~COperatorResultStore();

    Void  SetSpectrumName( const char* name );
    const std::string& GetSpectrumName() const;

    Void  StorePerTemplateResult( const CTemplate& t, const char* name, const COperatorResult& result );
    Void  StoreGlobalResult( const char* name, const COperatorResult& result );

    const COperatorResult*  GetPerTemplateResult( const CTemplate& t, const char* name ) const;
    TOperatorResultMap      GetPerTemplateResult( const char* name ) const;
    const COperatorResult*  GetGlobalResult( const char* name ) const;

    Void SaveAllResults( const char* dir ) const;

    Void PushScope( const char* name );
    Void PopScope();

    std::string GetCurrentScopeName() const;

protected:

    void StoreResult( TResultsMap& map, const char* name, const COperatorResult& result );

    void CreateResultStorage( std::fstream& stream, const boost::filesystem::path& path, const boost::filesystem::path& baseDir ) const;

    TPerTemplateResultsMap          m_PerTemplateResults;
    TResultsMap                     m_GlobalResults;

    std::string                     m_SpectrumName;

    TScopeStack                     m_ScopeStack;

};


}

#endif

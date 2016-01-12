#ifndef _REDSHIFT_PROCESSFLOW_DATASTORE_
#define _REDSHIFT_PROCESSFLOW_DATASTORE_

#include <epic/core/common/datatypes.h>

#include <epic/redshift/processflow/parameterstore.h>
#include <epic/redshift/processflow/resultstore.h>

#include <boost/filesystem.hpp>
#include <vector>
#include <ostream>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CDataStore
{

public:

    typedef std::vector<std::string>    TScopeStack;
    class CAutoScope {
        public:
            CAutoScope( CDataStore& store, const std::string& name );
            ~CAutoScope();
        private:
            CDataStore* m_Store;
    };

    CDataStore( COperatorResultStore& resultStore, CParameterStore& parameStore );
    virtual ~CDataStore();

    Void                PushScope( const std::string& name );
    Void                PopScope();

    const std::string&  GetSpectrumName() const;
    Void                SetSpectrumName( const std::string& name );

    std::string         GetCurrentScopeName() const;

    std::string         GetScope( const COperatorResult&  result) const;

    // Wrapper functions
    Bool                            GetScopedParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    Bool                            GetScopedParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    Bool                            GetScopedParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    Bool                            GetScopedParam( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    Bool                            GetScopedParam( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    Bool                            GetScopedParam( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    Bool                            GetScopedParam( const std::string& name, std::string& v, std::string defaultValue = "" ) const;

    Bool                            SetScopedParam( const std::string& name, const TFloat64List& v );
    Bool                            SetScopedParam( const std::string& name, const TInt64List& v );
    Bool                            SetScopedParam( const std::string& name, const TBoolList& v );
    Bool                            SetScopedParam( const std::string& name, Float64 v );
    Bool                            SetScopedParam( const std::string& name, Int64 v );
    Bool                            SetScopedParam( const std::string& name, Bool v );
    Bool                            SetScopedParam( const std::string& name, const std::string& v );

    Bool                            GetParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    Bool                            GetParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    Bool                            GetParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    Bool                            GetParam( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    Bool                            GetParam( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    Bool                            GetParam( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    Bool                            GetParam( const std::string& name, std::string& v, std::string defaultValue = "" ) const;

    Bool                            SetParam( const std::string& name, const TFloat64List& v );
    Bool                            SetParam( const std::string& name, const TInt64List& v );
    Bool                            SetParam( const std::string& name, const TBoolList& v );
    Bool                            SetParam( const std::string& name, Float64 v );
    Bool                            SetParam( const std::string& name, Int64 v );
    Bool                            SetParam( const std::string& name, Bool v );
    Bool                            SetParam( const std::string& name, const std::string& v );

    Void                            StoreScopedPerTemplateResult( const CTemplate& t, const std::string& name, std::shared_ptr<const COperatorResult>  result );
    Void                            StoreScopedGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult>  result );

    std::weak_ptr<const COperatorResult>          GetPerTemplateResult( const CTemplate& t, const std::string& name ) const;
    TOperatorResultMap              GetPerTemplateResult( const std::string& name ) const;
    std::weak_ptr<const COperatorResult>          GetGlobalResult( const std::string& name ) const;

    Void                            SaveRedshiftResultHeader( const boost::filesystem::path& dir );
    Void                            SaveRedshiftResult( const boost::filesystem::path& dir );
    Void                            SaveAllResults( const boost::filesystem::path& dir ) const;


protected:

    std::string                     GetScopedName( const std::string& name ) const;

    COperatorResultStore&            m_ResultStore;

    CParameterStore&                 m_ParameterStore;

    std::string                     m_SpectrumName;

    TScopeStack                     m_ScopeStack;


};


}

#endif

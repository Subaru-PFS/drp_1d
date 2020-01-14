#ifndef _REDSHIFT_PROCESSFLOW_DATASTORE_
#define _REDSHIFT_PROCESSFLOW_DATASTORE_

#include <RedshiftLibrary/common/datatypes.h>

#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/processflow/resultstore.h>

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

    void                PushScope( const std::string& name );
    void                PopScope();

    const std::string&  GetSpectrumName() const;
    void                SetSpectrumName( const std::string& name );
    const std::string&  GetProcessingID() const;
    void                SetProcessingID( const std::string& valStr );

    std::string         GetCurrentScopeName() const;

    std::string         GetScope( const COperatorResult&  result) const;

    // Wrapper functions
    void                            GetScopedParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    void                            GetScopedParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    void                            GetScopedParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    void                            GetScopedParam( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    void                            GetScopedParam( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    void                            GetScopedParam( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    void                            GetScopedParam( const std::string& name, std::string& v, std::string defaultValue = "" ) const;

    void                            SetScopedParam( const std::string& name, const TFloat64List& v );
    void                            SetScopedParam( const std::string& name, const TInt64List& v );
    void                            SetScopedParam( const std::string& name, const TBoolList& v );
    void                            SetScopedParam( const std::string& name, Float64 v );
    void                            SetScopedParam( const std::string& name, Int64 v );
    void                            SetScopedParam( const std::string& name, Bool v );
    void                            SetScopedParam( const std::string& name, const std::string& v );

    void                            GetParam( const std::string& name, TFloat64List& v, const TFloat64List& defaultValue = TFloat64List() ) const;
    void                            GetParam( const std::string& name, TInt64List& v, const TInt64List& defaultValue = TInt64List() ) const;
    void                            GetParam( const std::string& name, TBoolList& v, const TBoolList& defaultValue = TBoolList() ) const;
    void                            GetParam( const std::string& name, Float64& v, Float64 defaultValue  = 0 ) const;
    void                            GetParam( const std::string& name, Int64& v, Int64 defaultValue = 0 ) const;
    void                            GetParam( const std::string& name, Bool& v, Bool defaultValue = true ) const;
    void                            GetParam( const std::string& name, std::string& v, std::string defaultValue = "" ) const;

    void                            SetParam( const std::string& name, const TFloat64List& v );
    void                            SetParam( const std::string& name, const TInt64List& v );
    void                            SetParam( const std::string& name, const TBoolList& v );
    void                            SetParam( const std::string& name, Float64 v );
    void                            SetParam( const std::string& name, Int64 v );
    void                            SetParam( const std::string& name, Bool v );
    void                            SetParam( const std::string& name, const std::string& v );

    void                            StoreScopedPerTemplateResult( const CTemplate& t, const std::string& name, std::shared_ptr<const COperatorResult>  result );
    void                            StoreScopedGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult>  result );
    void                            ChangeScopedGlobalResult( const std::string& oldkey, const std::string& newkey );
    void                            DeleteScopedGlobalResult( const std::string& name );
    void                            StoreGlobalResult( const std::string& name, std::shared_ptr<const COperatorResult>  result );

    std::weak_ptr<const COperatorResult>          GetPerTemplateResult( const CTemplate& t, const std::string& name ) const;
    TOperatorResultMap              GetPerTemplateResult( const std::string& name ) const;
    std::weak_ptr<const COperatorResult>          GetGlobalResult( const std::string& name ) const;

    void                            SaveRedshiftResult( const boost::filesystem::path& dir );
    void                            SaveCandidatesResult( const boost::filesystem::path& dir );
    void                            SaveReliabilityResult( const boost::filesystem::path& dir );
    void                            SaveStellarResult( const boost::filesystem::path& dir );
    void                            SaveQsoResult( const boost::filesystem::path& dir );
    void                            SaveClassificationResult( const boost::filesystem::path& dir );
    void                            SaveAllResults(const boost::filesystem::path& dir , const std::string opt) const;


protected:

    std::string                     GetScopedName( const std::string& name ) const;

    COperatorResultStore&            m_ResultStore;

    CParameterStore&                 m_ParameterStore;

    std::string                     m_SpectrumName;
    std::string                     m_ProcessingID;

    TScopeStack                     m_ScopeStack;


};


}

#endif

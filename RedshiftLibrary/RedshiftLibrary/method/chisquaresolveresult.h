#ifndef _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARESOLVERESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CChisquareSolveResult : public COperatorResult
{

public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };

    CChisquareSolveResult();
    virtual ~CChisquareSolveResult();

    void SetType(const Int32 type);
    void SetScope(const std::string scope);
    void SetName(const std::string name);

    const Int32 GetType() const;
    const std::string GetScope() const;
    const std::string GetName() const;

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName, Float64 &amplitude, Float64 &amplitudeError, Float64 &dustCoeff, Int32 &meiksinIdx) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;
    Bool GetBestRedshiftFromPdf(const CDataStore& store, Float64& redshift, Float64& merit, Float64 &evidence) const;
    Int32 GetBestModel(const CDataStore& store, Float64 z, std::string& tplName) const;

    Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const;
    Bool GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates, Int32 n_candidates) const;

    Int32 m_type;
    std::string m_name;
    std::string m_scope;

    Int32 m_bestRedshiftMethod = 2; //best chi2, best proba

};


}

#endif

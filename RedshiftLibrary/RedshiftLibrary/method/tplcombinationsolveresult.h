#ifndef _REDSHIFT_OPERATOR_TPLCOMBINATIONSOLVERESULT_
#define _REDSHIFT_OPERATOR_TPLCOMBINATIONSOLVERESULT_

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
class CTplcombinationSolveResult : public COperatorResult
{

public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };

    CTplcombinationSolveResult();
    virtual ~CTplcombinationSolveResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName , Float64 &amplitude, Float64 &dustCoeff, Int32 &meiksinIdx) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;
    Bool GetBestRedshiftFromPdf(const CDataStore& store, Float64& redshift, Float64& merit, Float64 &evidence) const;
    Int32 GetBestModel(const CDataStore& store, Float64 z, std::string& tplName) const;

    Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const;
    Bool GetRedshiftCandidates( const CDataStore& store,  std::vector<Float64>& redshiftcandidates, Int32 n_candidates) const;

    Int32 m_type;

    Int32 m_bestRedshiftMethod = 2; //best chi2, best proba

};


}

#endif



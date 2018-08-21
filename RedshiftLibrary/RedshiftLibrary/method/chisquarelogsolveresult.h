#ifndef _REDSHIFT_OPERATOR_CHISQUARELOGSOLVERESULT_
#define _REDSHIFT_OPERATOR_CHISQUARELOGSOLVERESULT_

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
class CChisquareLogSolveResult : public COperatorResult
{

public:

    enum EType
    {
             nType_raw = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
             nType_all = 4,
    };

    CChisquareLogSolveResult();
    virtual ~CChisquareLogSolveResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit, std::string& tplName , Float64& amplitude, Float64& dustCoeff, Int32& meiksinIdx) const;
    Bool GetBestRedshiftPerTemplateString( const CDataStore& store, std::string& output ) const;
    Bool GetBestRedshiftFromPdf(const CDataStore& store, Float64& redshift, Float64& merit, Float64 &evidence) const;
    Int32 GetBestModel(const CDataStore& store, Float64 z, std::string& tplName) const;

    Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const;

    Int32 m_type;

    Int32 m_bestRedshiftMethod = 2; //0=bestChi2, 2=MargPDF
};


}

#endif



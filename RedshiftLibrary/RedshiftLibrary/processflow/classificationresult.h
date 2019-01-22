#ifndef _REDSHIFT_OPERATOR_CLASSIFICATIONRESULT_
#define _REDSHIFT_OPERATOR_CLASSIFICATIONRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/processflow/result.h>

#include <vector>
#include <ostream>
#include <map>

namespace NSEpic
{

class CDataStore;

/**
 * \ingroup Redshift
 */
class CClassificationResult : public COperatorResult
{

public:

    CClassificationResult();
    ~CClassificationResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }

    void SetTypeLabel( std::string lbl );
    void SetG(Float64 evidence);
    void SetS(Float64 evidence);
    void SetQ(Float64 evidence);


protected:

    std::string m_TypeLabel="-1";
    Float64 m_evidence_galaxy=-1.0;
    Float64 m_evidence_star=-1.0;
    Float64 m_evidence_qso=-1.0;
};

}

#endif

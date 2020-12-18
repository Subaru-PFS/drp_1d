#ifndef _REDSHIFT_PROCESSFLOW_CLASSIFICATIONRESULT_
#define _REDSHIFT_PROCESSFLOW_CLASSIFICATIONRESULT_

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
    void SetG(Float64 evidence, Float64 prob);
    void SetS(Float64 evidence, Float64 prob);
    void SetQ(Float64 evidence, Float64 prob);

  void getData(const std::string& name, std::string& v) const;
  void getData(const std::string& name, Float64& v) const;

protected:

    std::string m_TypeLabel="-1";
    Float64 m_evidence_galaxy=-1.0;
    Float64 m_evidence_star=-1.0;
    Float64 m_evidence_qso=-1.0;
    Float64 m_prob_galaxy=-1.0;
    Float64 m_prob_star=-1.0;
    Float64 m_prob_qso=-1.0;
};

}

#endif

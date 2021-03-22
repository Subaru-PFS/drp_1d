#ifndef _REDSHIFT_PROCESSFLOW_CLASSIFICATIONRESULT_
#define _REDSHIFT_PROCESSFLOW_CLASSIFICATIONRESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/solveresult.h>

#include <vector>
#include <ostream>
#include <map>

namespace NSEpic
{

class CDataStore;

/**
 * \ingroup Redshift
 */
class CClassificationResult : public CSolveResult
{

public:

    CClassificationResult();

    void Save( std::ostream& stream ) const;
    void SaveLine(std::ostream& stream ) const;

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

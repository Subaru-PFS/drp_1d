#ifndef _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE2_
#define _REDSHIFT_OPERATOR_LINEMATCHINGSOLVE2_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/linematching2solveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class COperatorLineMatching2Solve
{

public:

    COperatorLineMatching2Solve();
    ~COperatorLineMatching2Solve();

    std::shared_ptr<const CLineMatching2SolveResult> Compute( CDataStore& resultStore, 
							      const CSpectrum& spc, 
							      const TFloat64Range& lambdaRange, 
							      const TFloat64Range& redshiftsRange, 
							      Float64 redshiftStep, 
							      const CRayCatalog &restRayCatalog);
    const std::string GetDescription();

private:

    // Peak Detection
    Float64 m_winsize;
    Float64 m_minsize;
    Float64 m_maxsize;
    Float64 m_detectioncut;
    Float64 m_detectionnoiseoffset;
    Float64 m_cut;
    Float64 m_strongcut;
    Float64 m_enlargeRate;

    // Line Matching
    Int64 m_minMatchNum;
    Float64 m_tol;
    Bool m_dynamicCut;

};

}

#endif

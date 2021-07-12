#ifndef _REDSHIFT_METHOD_LINEMATCHINGSOLVE_
#define _REDSHIFT_METHOD_LINEMATCHINGSOLVE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/linematchingsolveresult.h"
#include "RedshiftLibrary/spectrum/template/template.h"

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 * \class CMethodLineMatchingSolve
 * \brief Solver method based on matching peaks to the lines catalogue.
 */
class CMethodLineMatchingSolve
{

public:

    CMethodLineMatchingSolve();
    ~CMethodLineMatchingSolve();

    std::shared_ptr<CLineMatchingSolveResult> Compute(CDataStore& resultStore,
                                                       const CSpectrum& spc, 
                                                       const TFloat64Range& lambdaRange, 
                                                       const TFloat64Range& redshiftsRange, 
                                                       Float64 redshiftStep, 
                                                       const CRayCatalog& restRayCatalog);

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
    Bool m_disablegaussianfitqualitycheck;
    Bool m_dynamicLinematching;
    Int64 m_minMatchNum;
    Float64 m_tol;

    // Log
    Bool m_bypassDebug; // If True, debug messages are suppressed even if the --verbose flag is passed.
};

}

#endif

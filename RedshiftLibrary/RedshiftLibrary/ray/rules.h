#ifndef _REDSHIFT_RAY_RULES_
#define _REDSHIFT_RAY_RULES_

#include <RedshiftLibrary/common/datatypes.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Rules
 */
class CRules
{

public:


    CRules(CSpectrum& spc, CRayCatalog& detectedCatalog, CRayCatalog& restCatalog, TFloat64Range &lambdaRange, Float64 winsize );
    ~CRules();


    Int32 check(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet);
private:

    Bool checkRule01(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet);
    Bool checkRule02(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet);
    Bool checkRule03(Float64 z, CRayMatchingResult::TSolutionSet& matchingSolutionSet);

    Float64 getRestRayLambda(std::string nametag);

    CSpectrum    m_spc;
    CRayCatalog    m_DetectedCatalog;
    CRayCatalog   m_RestCatalog;
    TFloat64Range m_lambdaRange;
    Float64 m_winsize;
};


}

#endif


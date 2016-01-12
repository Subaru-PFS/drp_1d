#ifndef _REDSHIFT_METHOD_CHISQUARESOLVE_
#define _REDSHIFT_METHOD_CHISQUARESOLVE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/chisquaresolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CMethodChisquareSolve
{

 public:

    enum EType
    {
             nType_full = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
    };

    CMethodChisquareSolve();
    ~CMethodChisquareSolve();
    std::shared_ptr<const CChisquareSolveResult> Compute(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold  , std::string opt_interp="lin");



private:

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold , Int32 spctype=nType_full, std::string opt_interp="lin");
};


}

#endif // _REDSHIFT_METHOD_CHISQUARESOLVE_


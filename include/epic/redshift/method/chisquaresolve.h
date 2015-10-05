#ifndef _REDSHIFT_METHOD_CHISQUARESOLVE_
#define _REDSHIFT_METHOD_CHISQUARESOLVE_


#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/chisquaresolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

class CMethodChisquareSolve : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CMethodChisquareSolve )

    public:

    enum EType
    {
             nType_full = 1,
             nType_continuumOnly = 2,
             nType_noContinuum = 3,
    };

    CMethodChisquareSolve();
    ~CMethodChisquareSolve();
    const CChisquareSolveResult *Compute( CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList,
                                        const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold  );



private:

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont,
                                   const TFloat64Range& lambdaRange, const TFloat64List& redshifts, Float64 overlapThreshold , Int32 spctype=nType_full);
};


}

#endif // _REDSHIFT_METHOD_CHISQUARESOLVE_


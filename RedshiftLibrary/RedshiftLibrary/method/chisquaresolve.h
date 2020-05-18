#ifndef _REDSHIFT_METHOD_CHISQUARESOLVE_
#define _REDSHIFT_METHOD_CHISQUARESOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/chisquaresolveresult.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>

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
             nType_noContinuum = 3
    };

    CMethodChisquareSolve( std::string calibrationPath );
    ~CMethodChisquareSolve();

    std::shared_ptr<const CChisquareSolveResult> Compute(CDataStore& resultStore,
                                                         const CSpectrum& spc,
                                                         const CTemplateCatalog& tplCatalog,
                                                         const TStringList& tplCategoryList,
                                                         const TFloat64Range& lambdaRange,
                                                         const TFloat64List& redshifts,
                                                         Float64 overlapThreshold,
                                                         const Float64 radius,
                                                         std::string opt_interp="lin");


private:

    Bool Solve(CDataStore& resultStore,
               const CSpectrum& spc,
               const CTemplate& tpl,
               const TFloat64Range& lambdaRange,
               const TFloat64List& redshifts,
               Float64 overlapThreshold,
               Int32 spctype=nType_full,
               std::string opt_interp="lin");

    std::string m_calibrationPath;
    Float64 m_radius;

};


}

#endif // _REDSHIFT_METHOD_CHISQUARESOLVE_

#ifndef _REDSHIFT_OPERATOR_LINEMODELTPLSHAPESOLVE_
#define _REDSHIFT_OPERATOR_LINEMODELTPLSHAPESOLVE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/method/linemodeltplshapesolveresult.h>
#include <epic/redshift/spectrum/template/template.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CLineModelTplshapeSolve
{

public:

    CLineModelTplshapeSolve();
    ~CLineModelTplshapeSolve();

    const std::string GetDescription();
    Bool PopulateParameters( CDataStore& dataStore );
    Bool LoadVelocities( const char* filePath, Float64& elv, Float64& alv );

    std::shared_ptr<const CLineModelTplshapeSolveResult> Compute( CDataStore& resultStore,
                                                                  const CSpectrum& spc,
                                                                  const CSpectrum& spcWithoutCont,
                                                                  const CTemplateCatalog &tplCatalog,
                                                                  const TStringList &tplCategoryList,
                                                                  const CRayCatalog& restraycatalog,
                                                                  const TFloat64Range& lambdaRange,
                                                                  const TFloat64List& redshifts );

    Bool Solve(CDataStore& resultStore,
               const CSpectrum& spc,
               const CSpectrum& spcCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList,
               const CTemplate& tpl,
               const CRayCatalog& lineCatalog,
               const TFloat64Range& lambdaRange,
               const TFloat64List& redshifts);

private:

    std::string m_opt_linetypefilter;
    std::string m_opt_lineforcefilter;
    std::string m_opt_fittingmethod;
    std::string m_opt_continuumcomponent;
    std::string m_opt_lineWidthType;
    Float64 m_opt_resolution;
    Float64 m_opt_velocity_emission;
    Float64 m_opt_velocity_absorption;
    std::string m_opt_velocityfit;
    std::string m_opt_continuumreest;
    std::string m_opt_rules;
    Float64 m_opt_extremacount;


};


}

#endif

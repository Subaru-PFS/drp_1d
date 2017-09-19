#ifndef _REDSHIFT_OPERATOR_LINEMODELMATCHINGSOLVE_
#define _REDSHIFT_OPERATOR_LINEMODELMATCHINGSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/linemodel.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CLineModelSolve
{

public:

    CLineModelSolve(std::string calibrationPath="");
    ~CLineModelSolve();

    const std::string GetDescription();
    Bool PopulateParameters( CDataStore& dataStore );


    std::shared_ptr<CLineModelSolveResult> Compute(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog& restraycatalog,
                                           const TFloat64Range& lambdaRange, const TFloat64List& redshifts );

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog& restraycatalog,
                                 const TFloat64Range& lambdaRange, const TFloat64List& redshifts);

private:

    Int32 CombinePDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result, std::string opt_rigidity, std::string opt_combine);
    Int32 SaveContinuumPDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result);


    std::string m_opt_linetypefilter;
    std::string m_opt_lineforcefilter;
    std::string m_opt_fittingmethod;
    std::string m_opt_continuumcomponent;
    std::string m_opt_rigidity;
    std::string m_opt_lineWidthType;
    Float64 m_opt_resolution;
    Float64 m_opt_velocity_emission;
    Float64 m_opt_velocity_absorption;
    std::string m_opt_velocityfit;
    Float64 m_opt_velocity_fit_min;
    Float64 m_opt_velocity_fit_max;
    std::string m_opt_continuumreest;
    std::string m_opt_rules;
    Float64 m_opt_extremacount;
    Float64 m_opt_twosteplargegridstep;
    std::string m_opt_combinePdf;
    std::string m_opt_pdfcombination;
    std::string m_opt_saveintermediateresults;
    bool m_opt_enableSaveChisquareTplshapeResults;

    std::string m_calibrationPath;
};


}

#endif

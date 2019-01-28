#ifndef _REDSHIFT_OPERATOR_LINEMODELMATCHINGSOLVE_
#define _REDSHIFT_OPERATOR_LINEMODELMATCHINGSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/linemodel.h>

#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/operator/pdfLogresult.h>

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
                                           const TFloat64Range& lambdaRange, const TFloat64List& redshifts , const string outputPdfRelDir);

    Bool Solve(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog& restraycatalog,
                                 const TFloat64Range& lambdaRange, const TFloat64List& redshifts);

private:

    Int32 CombinePDF(std::shared_ptr<const CLineModelResult> result,
                     std::string opt_rigidity,
                     std::string opt_combine,
                     Float64 opt_stronglinesprior,
                     Float64 opt_euclidNHaEmittersPriorStrength,
                     std::shared_ptr<CPdfMargZLogResult> postmargZResult,
                     std::shared_ptr<CPdfLogResult> zPrior);
    Int32 SaveContinuumPDF(CDataStore &store, std::shared_ptr<const CLineModelResult> result);


    std::string m_opt_linetypefilter;
    std::string m_opt_lineforcefilter;
    std::string m_opt_fittingmethod;
    std::string m_opt_secondpasslcfittingmethod;
    std::string m_opt_continuumcomponent;
    std::string m_opt_skipsecondpass="no";

    std::string m_opt_tplfit_dustfit="no";
    std::string m_opt_tplfit_igmfit="no";
    Float64 m_opt_continuumfitcount;
    std::string m_opt_tplfit_ignoreLinesSupport="no";

    std::string m_opt_rigidity;
    std::string m_opt_lineWidthType;
    Float64 m_opt_resolution;
    Float64 m_opt_velocity_emission;
    Float64 m_opt_velocity_absorption;
    std::string m_opt_velocityfit;
    Float64 m_opt_em_velocity_fit_min;
    Float64 m_opt_em_velocity_fit_max;
    Float64 m_opt_em_velocity_fit_step;
    Float64 m_opt_abs_velocity_fit_min;
    Float64 m_opt_abs_velocity_fit_max;
    Float64 m_opt_abs_velocity_fit_step;
    std::string m_opt_continuumreest;
    std::string m_opt_rules;

    //options for rigidity=tplshape
    std::string m_opt_tplratio_reldirpath="";
    std::string m_opt_tplratio_ismfit="no";
    std::string m_opt_offsets_reldirpath="";

    Float64 m_opt_extremacount;
    Float64 m_opt_candidatesLogprobaCutThreshold;
    Float64 m_opt_firstpass_largegridstep;
    std::string m_opt_firstpass_largegridsampling;
    std::string m_opt_firstpass_tplratio_ismfit;
    std::string m_opt_firstpass_disablemultiplecontinuumfit;
    std::string m_opt_firstpass_fittingmethod;

    std::string m_opt_pdfcombination;
    Float64 m_opt_stronglinesprior;
    Float64 m_opt_euclidNHaEmittersPriorStrength;
    std::string m_opt_saveintermediateresults;
    bool m_opt_enableSaveChisquareTplshapeResults;

    std::string m_calibrationPath;
};


}

#endif

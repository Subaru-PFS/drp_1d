#ifndef _REDSHIFT_METHOD_LINEMODELSOLVE_
#define _REDSHIFT_METHOD_LINEMODELSOLVE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/method/linemodelsolveresult.h>
#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/linemodel.h>
#include <RedshiftLibrary/operator/pdfz.h>
#include <RedshiftLibrary/operator/pdfMargZLogResult.h>
#include <RedshiftLibrary/operator/pdfLogresult.h>
#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/processflow/resultstore.h>

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;
class CDataStore;

/**
 * \ingroup Redshift
 */
class CLineModelSolve: public CSolve
{

public:

  CLineModelSolve(TScopeStack &scope,std::string objectType,std::string calibrationPath="");


    Bool PopulateParameters( std::shared_ptr<CParameterStore> parameterStore );

    std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                          std::shared_ptr<COperatorResultStore> resultStore,
                                          TScopeStack &scopeCDataStore);
  /*                                               const CSpectrum& spc,
                                                   const CTemplateCatalog& tplCatalog,
                                                   const TStringList& tplCategoryList,
                                                   const CRayCatalog& restraycatalog,
                                                   const TFloat64Range& lambdaRange,
                                                   const TFloat64List& redshifts,
                                                   const string outputPdfRelDir,
                                                   const Float64 radius);
  */

  Bool Solve(std::shared_ptr<COperatorResultStore> resultStore,
               const CSpectrum& spc,
               const CTemplateCatalog& tplCatalog,
               const TStringList& tplCategoryList,
               const CRayCatalog& restraycatalog,
               const TFloat64Range& lambdaRange,
               const TFloat64List& redshifts);

    ChisquareArray BuildChisquareArray(std::shared_ptr<const CLineModelResult> result,
                                        std::string opt_rigidity,
                                        std::string opt_combine,
                                        Float64 opt_stronglinesprior,
                                        Float64 opt_hapriorstrength,
                                        Float64 opt_euclidNHaEmittersPriorStrength,
                                        Float64 opt_modelPriorZStrength) const;

    //Int32 SaveContinuumPDF(CDataStore& store, std::shared_ptr<const CLineModelResult> result);

    void storeExtremaResults( std::shared_ptr<COperatorResultStore> dataStore,
                              std::shared_ptr<const CLineModelExtremaResult> ExtremaResult) const;

    void StoreChisquareTplShapeResults(std::shared_ptr<COperatorResultStore>  dataStore, std::shared_ptr<const CLineModelResult> result) const;


    COperatorLineModel m_linemodel;

    std::string m_opt_linetypefilter;
    std::string m_opt_lineforcefilter;
    std::string m_opt_enableLSF;
    std::string m_opt_fittingmethod;
    std::string m_opt_secondpasslcfittingmethod;
    std::string m_opt_continuumcomponent;
    std::string m_opt_skipsecondpass="no";
    std::string m_opt_secondpass_continuumfit="fromfirstpass";

    std::string m_opt_tplfit_method="templatefittinglog";
    std::string m_opt_tplfit_method_secondpass="templatefittinglog";
    std::string m_opt_tplfit_dustfit="no";
    std::string m_opt_tplfit_igmfit="no";
    Float64 m_opt_continuumfitcount;
    Float64 m_opt_tplfit_continuumprior_betaA=1.0;
    Float64 m_opt_tplfit_continuumprior_betaTE=1.0;
    Float64 m_opt_tplfit_continuumprior_betaZ=1.0;
    Float64 m_opt_continuum_neg_amp_threshold=-INFINITY; // no thresholding
    std::string m_opt_tplfit_continuumprior_dirpath="";
    std::string m_opt_tplfit_ignoreLinesSupport="no";

    std::string m_opt_rigidity;
    std::string m_opt_lineWidthType;
    Float64 m_opt_nsigmasupport;
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
    std::string m_opt_enableImproveBalmerFit;

    std::string m_opt_lya_forcefit;
    std::string m_opt_lya_forcedisablefit;
    Float64 m_opt_lya_fit_asym_min;
    Float64 m_opt_lya_fit_asym_max;
    Float64 m_opt_lya_fit_asym_step;
    Float64 m_opt_lya_fit_width_min;
    Float64 m_opt_lya_fit_width_max;
    Float64 m_opt_lya_fit_width_step;
    Float64 m_opt_lya_fit_delta_min;
    Float64 m_opt_lya_fit_delta_max;
    Float64 m_opt_lya_fit_delta_step;


    //options for rigidity=tplshape
    std::string m_opt_tplratio_reldirpath="";
    std::string m_opt_tplratio_ismfit="no";
    Float64 m_opt_tplratio_prior_betaA=1.0;
    Float64 m_opt_tplratio_prior_betaTE=1.0;
    Float64 m_opt_tplratio_prior_betaZ=1.0;
    std::string m_opt_tplratio_prior_dirpath="";
    std::string m_opt_offsets_reldirpath="";

    Int64 m_opt_extremacount;
    Int64 m_opt_extremacountB;

    Float64 m_opt_candidatesLogprobaCutThreshold;
    Float64 m_opt_firstpass_largegridstep;
    std::string m_opt_firstpass_largegridsampling;
    std::string m_opt_firstpass_tplratio_ismfit;
    std::string m_opt_firstpass_disablemultiplecontinuumfit;
    std::string m_opt_firstpass_fittingmethod;

    std::string m_opt_pdfcombination;
    std::string m_opt_pdf_margAmpCorrection="no";
    Float64 m_opt_stronglinesprior;
    Float64 m_opt_haPrior;
    Float64 m_opt_euclidNHaEmittersPriorStrength;
    Float64 m_opt_modelZPriorStrength;
    std::string m_opt_saveintermediateresults;
    bool m_opt_enableSaveChisquareTplshapeResults;
    Float64 m_opt_secondpass_halfwindowsize;

    std::string m_calibrationPath;
    std::string m_outputPdfRelDir;
    Float64 m_redshiftSeparation;

};


}

#endif

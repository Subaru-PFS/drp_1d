// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_METHOD_LINEMODELSOLVE_
#define _REDSHIFT_METHOD_LINEMODELSOLVE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/method/linemodelsolveresult.h"
#include "RedshiftLibrary/method/solve.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/operator/linemodel.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/operator/pdfMargZLogResult.h"
#include "RedshiftLibrary/operator/pdfLogresult.h"
#include "RedshiftLibrary/processflow/inputcontext.h"
#include "RedshiftLibrary/processflow/resultstore.h"

namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;

/**
 * \ingroup Redshift
 */
class CLineModelSolve: public CSolve
{
public:

    CLineModelSolve(TScopeStack &scope,std::string objectType,std::string calibrationPath="");

    bool PopulateParameters( std::shared_ptr<const CParameterStore> parameterStore );

    std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                          std::shared_ptr<COperatorResultStore> resultStore,
                                          TScopeStack &scope) override;

    bool Solve(std::shared_ptr<COperatorResultStore> resultStore,
               const CSpectrum& spc,
               const CSpectrum& rebinnedSpc,
               const CTemplateCatalog& tplCatalog,
               const TStringList& tplCategoryList,
               const CRayCatalog::TRayVector& restraycatalog,
               const CRayCatalogsTplShape& tplRatioCatalog,
               const TFloat64Range& lambdaRange,
               const TFloat64List& redshifts,
               const std::shared_ptr<const CPhotBandCatalog> &photBandCat,
               const Float64 photo_weight);

private:

    ChisquareArray BuildChisquareArray(std::shared_ptr<const CLineModelResult> result,
                                        std::string opt_rigidity,
                                        std::string opt_combine,
                                        Float64 opt_stronglinesprior,
                                        Float64 opt_hapriorstrength,
                                        Float64 opt_euclidNHaEmittersPriorStrength) const;

   
      void storeExtremaResults( std::shared_ptr<COperatorResultStore> dataStore,
                             std::shared_ptr<const LineModelExtremaResult> ExtremaResult) const;

    void StoreChisquareTplShapeResults(std::shared_ptr<COperatorResultStore>  dataStore, std::shared_ptr<const CLineModelResult> result) const;


    COperatorLineModel m_linemodel;

    std::string m_opt_linetypefilter;
    std::string m_opt_lineforcefilter;
    std::string m_opt_fittingmethod;
    std::string m_opt_secondpasslcfittingmethod;
    std::string m_opt_continuumcomponent;
    bool m_opt_skipsecondpass=false;
    std::string m_opt_secondpass_continuumfit="fromfirstpass";

    bool m_opt_tplfit_fftprocessing=true;//default to using fft
    bool m_opt_tplfit_fftprocessing_secondpass=true;
    bool m_opt_tplfit_use_photometry=false;
    Float64 m_opt_tplfit_photo_weight = 1.0;
    bool m_opt_tplfit_dustfit=false;
    bool m_opt_tplfit_igmfit=false;
  Float64 m_opt_continuumfitcount; //TODO is this really a double and not an integer ?
    Float64 m_opt_tplfit_continuumprior_betaA=1.0;
    Float64 m_opt_tplfit_continuumprior_betaTE=1.0;
    Float64 m_opt_tplfit_continuumprior_betaZ=1.0;
    Float64 m_opt_continuum_neg_amp_threshold=-INFINITY; // no thresholding
    std::string m_opt_tplfit_continuumprior_dirpath="";
    bool m_opt_tplfit_ignoreLinesSupport=false;

    std::string m_opt_rigidity;
    std::string m_opt_lineWidthType;
    Float64 m_opt_nsigmasupport;
    Float64 m_opt_velocity_emission;
    Float64 m_opt_velocity_absorption;
    bool m_opt_velocityfit;
    Float64 m_opt_em_velocity_fit_min;
    Float64 m_opt_em_velocity_fit_max;
    Float64 m_opt_em_velocity_fit_step;
    Float64 m_opt_abs_velocity_fit_min;
    Float64 m_opt_abs_velocity_fit_max;
    Float64 m_opt_abs_velocity_fit_step;
    std::string m_opt_continuumreest;
    std::string m_opt_rules;
    bool m_opt_enableImproveBalmerFit;

    bool m_opt_lya_forcefit;
    bool m_opt_lya_forcedisablefit;
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
    bool m_opt_tplratio_ismfit=false;
    Float64 m_opt_tplratio_prior_betaA=1.0;
    Float64 m_opt_tplratio_prior_betaTE=1.0;
    Float64 m_opt_tplratio_prior_betaZ=1.0;
    std::string m_opt_tplratio_prior_dirpath="";
    std::string m_opt_offsets_reldirpath="";

    Int64 m_opt_extremacount;
    Int64 m_opt_extremacountB;

    Float64 m_opt_candidatesLogprobaCutThreshold;
    UInt32 m_opt_firstpass_largegridstepRatio;
    std::string m_opt_firstpass_largegridsampling;
    bool m_opt_firstpass_tplratio_ismfit;
    bool m_opt_firstpass_disablemultiplecontinuumfit;
    std::string m_opt_firstpass_fittingmethod;

    std::string m_opt_pdfcombination;
    bool m_opt_pdf_margAmpCorrection=false;
    Float64 m_opt_stronglinesprior;
    Float64 m_opt_haPrior;
    Float64 m_opt_euclidNHaEmittersPriorStrength;
    std::string m_opt_saveintermediateresults;
    Float64 m_opt_secondpass_halfwindowsize;

    std::string m_calibrationPath;
    std::string m_outputPdfRelDir;
    Float64 m_redshiftSeparation;

};


}

#endif

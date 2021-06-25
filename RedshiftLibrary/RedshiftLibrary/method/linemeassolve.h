#ifndef _LINEMEAS_SOLVE_H
#define _LINEMEAS_SOLVE_H

#include <RedshiftLibrary/common/datatypes.h>

#include <RedshiftLibrary/method/solve.h>
#include <RedshiftLibrary/method/linemeassolveresult.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/operator/linemodel.h>
#include <RedshiftLibrary/processflow/context.h>


namespace NSEpic
{

class CSpectrum;
class CTemplateCatalog;


/**
 * \ingroup Redshift
 */
  class CLineMeasSolve : public CSolve
{

public:

  CLineMeasSolve(TScopeStack &scope,std::string objectType,std::string calibrationPath="");
  ~CLineMeasSolve();

  void solve();
  void Init();
  std::shared_ptr<CSolveResult> compute(std::shared_ptr<const CInputContext> inputContext,
                                        std::shared_ptr<COperatorResultStore> resultStore,
                                        TScopeStack &scope);

  void GetRedshiftSampling(std::shared_ptr<const CInputContext> inputContext, TFloat64Range& redshiftRange, Float64& redshiftStep); 
protected:
  COperatorLineModel m_linemodel;

  std::string m_calibrationPath;

  /*
  TRedshiftList m_redshiftRange;
  Float64 m_zref;

  std::string m_opt_fittingmethod;
  std::string m_opt_rigidity;

  std::string m_opt_lineWidthType;
  Float64 m_opt_velocity_emission;
  Float64 m_opt_velocity_absorption;
  std::string m_opt_velocityfit;
  Float64 m_opt_em_velocity_fit_min;
  Float64 m_opt_em_velocity_fit_max;
  Float64 m_opt_em_velocity_fit_step;
  Float64 m_opt_abs_velocity_fit_min;
  Float64 m_opt_abs_velocity_fit_max;
  Float64 m_opt_abs_velocity_fit_step;
  Float64 m_opt_manvelfit_dz_min;
  Float64 m_opt_manvelfit_dz_max;
  Float64 m_opt_manvelfit_dz_step;
  std::string m_opt_continuumreest;
  std::string m_opt_rules;
  //  std::string m_opt_enableImproveBalmerFit;
  std::string m_opt_secondpass_continuumfit;
  */
};

}
#endif

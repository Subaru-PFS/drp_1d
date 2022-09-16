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
#ifndef _LINEMEAS_SOLVE_H
#define _LINEMEAS_SOLVE_H

#include <RedshiftLibrary/common/datatypes.h>

#include <RedshiftLibrary/method/linemeassolveresult.h>
#include <RedshiftLibrary/method/objectSolve.h>
#include <RedshiftLibrary/operator/linemodel.h>
#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/spectrum/template/template.h>

namespace NSEpic {

class CSpectrum;
class CTemplateCatalog;

/**
 * \ingroup Redshift
 */
class CLineMeasSolve : public CObjectSolve {

public:
  CLineMeasSolve(TScopeStack &scope, std::string objectType);

  void Init();
  std::shared_ptr<CSolveResult>
  compute(std::shared_ptr<const CInputContext> inputContext,
          std::shared_ptr<COperatorResultStore> resultStore,
          TScopeStack &scope);

  void GetRedshiftSampling(std::shared_ptr<const CInputContext> inputContext,
                           TFloat64Range &redshiftRange, Float64 &redshiftStep);

protected:
  COperatorLineModel m_linemodel;

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

} // namespace NSEpic
#endif

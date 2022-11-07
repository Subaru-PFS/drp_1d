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
#include "RedshiftLibrary/operator/logZPdfResult.h"
#include "RedshiftLibrary/common/vectorOperations.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <string>

#include "RedshiftLibrary/log/log.h"
#include <boost/algorithm/string/predicate.hpp>
#include <gsl/gsl_interp.h>
#include <iomanip>
#include <iostream>
using namespace std;
using namespace NSEpic;

CLogZPdfResult::CLogZPdfResult() { this->m_type = "CLogZPdfResult"; }
CLogZPdfResult::CLogZPdfResult(const TZGridListParams &zparamList,
                               const TFloat64List &valproba_mixed)
    : valProbaLog(valproba_mixed), Redshifts(valproba_mixed.size(), -DBL_MAX) {
  this->m_type = "CLogZPdfResult";
  setZGridParams(zparamList);
}
CLogZPdfResult::CLogZPdfResult(const TFloat64List &redshifts,
                               const TZGridListParams &zparamList)
    : Redshifts(redshifts), valProbaLog(redshifts.size(), -DBL_MAX) {
  this->m_type = "CLogZPdfResult";
  setZGridParams(zparamList); // to decide if we replace this
}

void CLogZPdfResult::setZGridParams(const TZGridListParams &paramList) {
  Int32 s = paramList.size();
  zcenter.resize(s, NAN);
  zmin.resize(s, NAN);
  zmax.resize(s, NAN);
  zstep.resize(s, NAN);

  std::transform(paramList.begin(), paramList.end(), zcenter.begin(),
                 [](const ZGridParameters &params) { return params.zcenter; });
  std::transform(paramList.begin(), paramList.end(), zmin.begin(),
                 [](const ZGridParameters &params) { return params.zmin; });
  std::transform(paramList.begin(), paramList.end(), zmax.begin(),
                 [](const ZGridParameters &params) { return params.zmax; });
  std::transform(paramList.begin(), paramList.end(), zstep.begin(),
                 [](const ZGridParameters &params) { return params.zstep; });
}

// usefull in the case of skipsecondpass, TemplateFitting and TplCombination
bool CLogZPdfResult::isZGridCoherent(const TFloat64List &zstep) {
  if (!zstep.size())
    THROWG(INTERNAL_ERROR, "zstep is empty");
  if (zstep.size() == 1)
    return true;
  // TODO: check if Redshifts grid has a coherent step, using a different method
  /*auto it = std::find_if(zstep.cbegin(), zstep.cend(),
                         [&](Float64 c) { c != zstep[0]; });
  bool isCoherent = it == zstep.end();
  return isCoherent;*/
  for (auto a : zstep)
    if (a != zstep[0])
      return false;

  return true;
}

// common function to create either coarse or fine grid, default to coarse
const TFloat64List CLogZPdfResult::buildLogZPdfGrid(
    bool logsampling, const TZGridListParams &zparams, ZGridType zgrid_type) {

  Float64 zgridStep = zparams.front().zstep;
  TFloat64Range zrange(zparams.front().zmin, zparams.front().zmax);

  switch (zgrid_type) {
  case FINE: {
    TFloat64List zsteps(zparams.size());
    std::transform(zparams.begin(), zparams.end(), zsteps.begin(),
                   [](const ZGridParameters &params) { return params.zstep; });
    if (!CLogZPdfResult::isZGridCoherent(zsteps))
      zgridStep = zparams[1].zstep;
    break;
  }
  case MIXED:
    return CLogZPdfResult::buildLogMixedZPdfGrid(logsampling, zparams);
  case COARSE:
  default:
    break;
  }

  return logsampling ? zrange.SpreadOverLogZplusOne(zgridStep)
                     : zrange.SpreadOver(zgridStep);
}

// mixed pdf
const TFloat64List
CLogZPdfResult::buildLogMixedZPdfGrid(bool logsampling,
                                      const TZGridListParams &zparams) {
  TFloat64List mixedGrid =
      CLogZPdfResult::buildLogZPdfGrid(logsampling, zparams, COARSE);
  // insert fineGrid intervals
  for (Int32 i = 1; i < zparams.size(); i++) {
    TFloat64Range range(zparams[i].zmin, zparams[i].zmax);
    Int32 imin = -1;
    Int32 imax = -1;
    bool b = range.getClosedIntervalIndices(mixedGrid, imin, imax, false);
    Int32 ndup = imax - imin + 1;
    const TFloat64List &extendedZ =
        logsampling ? range.SpreadOverLogZplusOne(zparams[i].zstep)
                    : range.SpreadOver(zparams[i].zstep);
    insertWithDuplicates(mixedGrid, imin, extendedZ, ndup);
  }
  return mixedGrid;
}

const TPdf CLogZPdfResult::getLogZPdf_fine(bool logsampling,
                                           const TZGridListParams &zparams,
                                           const TFloat64List &valProbaLog) {
  TPdf res;
  TFloat64List zsteps(zparams.size());
  std::transform(zparams.begin(), zparams.end(), zsteps.begin(),
                 [](const ZGridParameters &params) { return params.zstep; });

  ZGridType zgtype = CLogZPdfResult::isZGridCoherent(zsteps) ? COARSE : MIXED;
  TFloat64List origin_zgrid =
      CLogZPdfResult::buildLogZPdfGrid(logsampling, zparams, zgtype);

  if (CLogZPdfResult::isZGridCoherent(zsteps) && zsteps.size() == 1)
    THROWG(INTERNAL_ERROR, "Cannot create fine grid out of coarse "
                           "parameters"); // return
                                          // TPdf(origin_zgrid,valProbaLog);

  res.zgrid = CLogZPdfResult::buildLogZPdfGrid(logsampling, zparams, FINE);

  interpolateLargeGridOnFineGrid(origin_zgrid, res.zgrid, valProbaLog,
                                 res.probaLog);
  return res;
}

void CLogZPdfResult::interpolateLargeGridOnFineGrid(
    const TFloat64List &originGrid, const TFloat64List &targetGrid,
    const TFloat64List &originValues, TFloat64List &outputValues) {
  Log.LogDetail("interp FROM large grid z0=%f to zEnd=%f (n=%d)",
                originGrid.front(), originGrid.back(), originGrid.size());
  Log.LogDetail("interp TO fine grid z0=%f to zEnd=%f (n=%d)",
                targetGrid.front(), targetGrid.back(), targetGrid.size());

  outputValues.resize(targetGrid.size(), NAN);
  // initialise and allocate the gsl objects
  // lin
  gsl_interp *interpolation =
      gsl_interp_alloc(gsl_interp_linear, originValues.size());
  gsl_interp_init(interpolation, &(originGrid.front()), &(originValues.front()),
                  originValues.size());
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

  for (Int32 j = 0; j < targetGrid.size(); j++) {
    Float64 Xrebin = targetGrid[j];
    if (Xrebin < originGrid.front() || Xrebin > originGrid.back()) {
      continue;
    }
    outputValues[j] =
        gsl_interp_eval(interpolation, &originGrid.front(),
                        &originValues.front(), Xrebin, accelerator); // lin
  }

  gsl_interp_free(interpolation);
  gsl_interp_accel_free(accelerator);
}

/*
TPdf CLogZPdfResult::getFineGridPdf(const TZGridListParams& pdf_params,
                                    const Float64* pdf_probaLog, Int32 n)
{
  TFloat64List ppl(n);
  for(unsigned int i=0;i<n;i++) ppl[i]=pdf_probaLog[i];
  CLogZPdfResult pdf_result(pdf_params,ppl);
  CLogZPdfResult fine = pdf_result.getLogZPdf_fine(true);
  return TPdf(fine.Redshifts,fine.valProbaLog);
}
*/

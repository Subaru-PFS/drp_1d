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
#include <algorithm>

#include <gsl/gsl_interp.h>

#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/common/vectorOperations.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/operator/logZPdfResult.h"
#include "RedshiftLibrary/operator/pdfz.h"

using namespace NSEpic;

CLogZPdfResult::CLogZPdfResult(const CZGridListParams &zparamList,
                               bool logsampling_, TFloat64List valProbaLog_)
    : COperatorResult("CLogZPdfResult"), valProbaLog(std::move(valProbaLog_)),
      logsampling(logsampling_) {
  setZGridParams(zparamList);
  if (valProbaLog.empty())
    valProbaLog.assign(redshifts.size(), -DBL_MAX);
  check_sizes();
}

void CLogZPdfResult::check_sizes() const {
  if (valProbaLog.size() != redshifts.size())
    THROWG(ErrorCode::INTERNAL_ERROR,
           Formatter() << "redshifts size (" << redshifts.size()
                       << ") is different from pdf size (" << valProbaLog.size()
                       << ")");
}

void CLogZPdfResult::setZGrid() {
  redshifts = getZGridParams().getZGrid(logsampling);
}

void CLogZPdfResult::setZGridParams(const CZGridListParams &paramList) {
  Int32 s = paramList.size();
  zcenter.resize(s, NAN);

  zmin.resize(s, NAN);
  zmax.resize(s, NAN);
  zstep.resize(s, NAN);

  std::transform(paramList.cbegin(), paramList.cend(), zcenter.begin(),
                 [](const CZGridParam &param) { return param.zcenter; });
  std::transform(paramList.cbegin(), paramList.cend(), zmin.begin(),
                 [](const CZGridParam &param) { return param.zmin; });
  std::transform(paramList.cbegin(), paramList.cend(), zmax.begin(),
                 [](const CZGridParam &param) { return param.zmax; });
  std::transform(paramList.cbegin(), paramList.cend(), zstep.begin(),
                 [](const CZGridParam &param) { return param.zstep; });
  setZGrid();
}

CZGridListParams CLogZPdfResult::getZGridParams() const {
  TZGridListParams zparams(zmin.size());
  for (int i = 0; i < ssize(zcenter); i++)
    zparams[i] =
        CZGridParam(TFloat64Range(zmin[i], zmax[i]), zstep[i], zcenter[i]);
  return CZGridListParams(std::move(zparams));
}

void CLogZPdfResult::interpolateOnGrid(TFloat64List targetGrid) {
  Log.LogDetail(Formatter()
                << "interp FROM initial grid z0=" << redshifts.front()
                << " to zEnd=" << redshifts.back() << " (n=" << redshifts.size()
                << ")");
  Log.LogDetail(Formatter() << "interp TO final grid z0=" << targetGrid.front()
                            << " to zEnd=" << targetGrid.back()
                            << " (n=" << targetGrid.size() << ")");

  TFloat64List outputValues(targetGrid.size(), NAN);

  // initialise and allocate the gsl objects
  // lin
  gsl_interp *interpolation =
      gsl_interp_alloc(gsl_interp_linear, valProbaLog.size());
  gsl_interp_init(interpolation, &(redshifts.front()), &(valProbaLog.front()),
                  valProbaLog.size());
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

  for (Int32 j = 0; j < ssize(targetGrid); j++) {
    Float64 Xrebin = targetGrid[j];
    if (Xrebin < redshifts.front() || Xrebin > redshifts.back()) {
      continue;
    }
    outputValues[j] =
        gsl_interp_eval(interpolation, &redshifts.front(), &valProbaLog.front(),
                        Xrebin, accelerator); // lin
  }

  gsl_interp_free(interpolation);
  gsl_interp_accel_free(accelerator);

  // update members
  redshifts = std::move(targetGrid);
  valProbaLog = std::move(outputValues);

  check_sizes();
}

void CLogZPdfResult::convertToRegular(bool fine) {

  if (isRegular())
    return;

  // set new zgrid params
  if (fine) {
    for (Int32 i = 2; i < ssize(zstep); ++i)
      if (zstep[1] != zstep[i])
        THROWG(ErrorCode::INTERNAL_ERROR,
               "cannot convert to regular with different steps "
               "in 2nd pass subgrids");
    zstep[0] = zstep[1];
  }
  zcenter.resize(1);
  zmin.resize(1);
  zmax.resize(1);
  zstep.resize(1);

  auto param = getZGridParams()[0];
  if (fine)
    param.zstep = zstep[1];
  auto target_redshifts = param.getZGrid(logsampling);

  interpolateOnGrid(target_redshifts);
}

void CLogZPdfResult::extrapolate_on_right_border(Float64 zend) {
  if (zend <= zmax[0])
    return;

  zmax[0] = zend;
  setZGrid();
  valProbaLog.resize(redshifts.size(), valProbaLog.back());
  check_sizes();
}

Float64 CLogZPdfResult::getSumTrapez() const {

  check_sizes();

  Float64 sum = 0.0;
  if (redshifts.size() < 2)
    return sum;

  // prepare LogEvidence
  Int32 sumMethod = 1;
  Float64 logSum = COperatorPdfz::logSumExpTrick(valProbaLog, redshifts);
  sum = exp(logSum);

  return sum;
}

void CLogZPdfResult::checkPdfSum() const {

  // check pdf sum=1
  Float64 sumTrapez = getSumTrapez();
  Log.LogDetail(
      Formatter()
      << "CLogZPdfResult::checkPdfSum: Pdfz normalization - sum trapz. = "
      << sumTrapez);
  if (abs(sumTrapez - 1.0) > 1e-1)
    THROWG(ErrorCode::PDF_NORMALIZATION_FAILED,
           Formatter() << "Pdfz normalization failed, "
                          "trapzesum = "
                       << sumTrapez);
}

void CLogZPdfResult::isPdfValid() const {

  if (redshifts.size() < 2)
    THROWG(ErrorCode::INTERNAL_ERROR, "PDF has size less than 2");

  // is it completely flat ?
  const auto minmax_it =
      std::minmax_element(valProbaLog.cbegin(), valProbaLog.cend());
  const Float64 minVal = *minmax_it.first;
  const Float64 maxVal = *minmax_it.second;

  if (minVal == maxVal)
    THROWG(ErrorCode::FLAT_ZPDF, "PDF is flat");

  // is pdf any value nan ?
  for (Int32 k = 0; k < ssize(valProbaLog); k++)
    if (valProbaLog[k] != valProbaLog[k])
      THROWG(ErrorCode::INTERNAL_ERROR, "PDF has nan or invalid values");

  // is sum equal to 1
  checkPdfSum();
}

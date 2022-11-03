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
bool CLogZPdfResult::isZGridCoherent() const {
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
const TFloat64List CLogZPdfResult::buildLogZPdfGrid(bool logsampling, bool fine,
                                                    Int32 idx) const {
  if (idx == -1 || idx > zmin.size() - 1)
    THROWG(INTERNAL_ERROR, "Index out of range");
  Float64 zgridStep = zstep.front();
  TFloat64Range zrange(zmin[idx], zmax[idx]);
  if (fine && !isZGridCoherent())
    zgridStep = zstep[1];
  return logsampling ? zrange.SpreadOverLogZplusOne(zgridStep)
                     : zrange.SpreadOver(zgridStep);
}

// mixed pdf
const TFloat64List
CLogZPdfResult::buildLogMixedZPdfGrid(bool logsampling) const {
  TFloat64List mixedGrid = buildLogZPdfGrid(logsampling, false);
  // insert fineGrid intervals
  for (Int32 i = 1; i < zmin.size(); i++) {
    TFloat64Range range(zmin[i], zmax[i]);
    Int32 imin = -1;
    Int32 imax = -1;
    bool b = range.getClosedIntervalIndices(mixedGrid, imin, imax, false);
    Int32 ndup = imax - imin + 1;
    const TFloat64List &extendedZ = buildLogZPdfGrid(logsampling, true, i);
    insertWithDuplicates(mixedGrid, imin, extendedZ, ndup);
  }
  return mixedGrid;
}

const CLogZPdfResult CLogZPdfResult::getLogZPdf_fine(bool logsampling) const {
  if (isZGridCoherent())
    return *this;
  CLogZPdfResult fineZPdf;
  fineZPdf.Redshifts = buildLogZPdfGrid(logsampling, true);
  fineZPdf.zmin.push_back(zmin.front());
  fineZPdf.zmax.push_back(zmax.front());
  fineZPdf.zstep.push_back(zstep[1]);
  fineZPdf.zcenter.push_back(NAN);

  interpolateLargeGridOnFineGrid(Redshifts, fineZPdf.Redshifts, valProbaLog,
                                 fineZPdf.valProbaLog);
  return fineZPdf;
}

void CLogZPdfResult::interpolateLargeGridOnFineGrid(
    const TFloat64List &coarseGrid, const TFloat64List &fineGrid,
    const TFloat64List &entityLargeGrid, TFloat64List &entityFineGrid) {
  Log.LogDetail("interp FROM large grid z0=%f to zEnd=%f (n=%d)",
                coarseGrid.front(), coarseGrid.back(), coarseGrid.size());
  Log.LogDetail("interp TO fine grid z0=%f to zEnd=%f (n=%d)", fineGrid.front(),
                fineGrid.back(), fineGrid.size());

  entityFineGrid.resize(fineGrid.size(), NAN);
  // initialise and allocate the gsl objects
  // lin
  gsl_interp *interpolation =
      gsl_interp_alloc(gsl_interp_linear, entityLargeGrid.size());
  gsl_interp_init(interpolation, &(coarseGrid.front()),
                  &(entityLargeGrid.front()), entityLargeGrid.size());
  gsl_interp_accel *accelerator = gsl_interp_accel_alloc();

  for (Int32 j = 0; j < fineGrid.size(); j++) {
    Float64 Xrebin = fineGrid[j];
    if (Xrebin < coarseGrid.front() || Xrebin > coarseGrid.back()) {
      continue;
    }
    entityFineGrid[j] =
        gsl_interp_eval(interpolation, &coarseGrid.front(),
                        &entityLargeGrid.front(), Xrebin, accelerator); // lin
  }

  gsl_interp_free(interpolation);
  gsl_interp_accel_free(accelerator);
}


TPdf CLogZPdfResult::getFineGridPdf(const TZGridListParams& pdf_params,
				    const TFloat64List& pdf_probaLog,
				    bool logsampling)
{
  CLogZPdfResult pdf_result(pdf_params,pdf_probaLog);
  CLogZPdfResult fine = pdf_result.getLogZPdf_fine(logsampling);
  return TPdf(fine.Redshifts,fine.valProbaLog);
}

TPdf CLogZPdfResult::getMixedGridPdf(const TZGridListParams& pdf_params,
				    const TFloat64List& pdf_probaLog,
				     bool logsampling)
{
  CLogZPdfResult pdf_result(pdf_params,pdf_probaLog);
  TFloat64List fine = pdf_result.buildLogMixedZPdfGrid(logsampling);
  return TPdf(fine,pdf_probaLog);
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

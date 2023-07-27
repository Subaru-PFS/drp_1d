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
#ifndef _REDSHIFT_LINE_CATALOGSTPLRATIO_
#define _REDSHIFT_LINE_CATALOGSTPLRATIO_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/lineRatioCatalog.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

#include <boost/format.hpp>

#include <string>
#include <vector>

namespace NSEpic {
class CLineModelElementList;

/**
 * \ingroup Redshift
 */
class CLineCatalogsTplRatio {

public:
  Float64 GetBestFit(const CLineMap &restLineList,
                     const TFloat64Map &fittedAmplitudes,
                     const TFloat64Map &fittedErrors,
                     TFloat64Map &amplitudesCorrected,
                     std::string &bestTplName) const;
  CLineMap GetRestLinesList(Int32 index) const;
  Int32 GetCatalogsCount() const { return m_lineRatioCatalogs.size(); }

  TFloat64List getCatalogsPriors() const;
  std::string GetCatalogName(Int32 idx) const {
    return m_lineRatioCatalogs.at(idx).getName();
  }
  Int32 GetIsmIndex(Int32 idx) const {
    return m_lineRatioCatalogs.at(idx).getIsmIndex();
  }

  bool GetCatalogVelocities(Int32 idx, Float64 &elv, Float64 &alv) const;
  std::vector<TFloat64Map> InitLineCorrespondingAmplitudes(
      Int32 enableISMCalzetti,
      const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
          &ismCorrectionCalzetti) const;
  const CLineRatioCatalog &GetCatalog(Int32 iCatalog) const {
    return m_lineRatioCatalogs[iCatalog];
  }

  void addLineRatioCatalog(const CLineRatioCatalog &lr_catalog) {
    m_lineRatioCatalogs.push_back(lr_catalog);
  }

private:
  Float64 computeFitValue(const TFloat64Map &ampsLM, const TFloat64Map &errLM,
                          const TFloat64Map &ampsTPL,
                          TFloat64Map &ampsCorrected) const;
  Float64 getFitForOneCatalog(const CLineMap &restLineList,
                              const TFloat64Map &fittedAmplitudes,
                              const TFloat64Map &fittedErrors,
                              const CLineRatioCatalog &catalog,
                              TFloat64Map &ampsCorrected) const;
  void logLineNominalAmp(
      bool enableISMCalzetti,
      const std::vector<TFloat64Map> &lineCatalogLinesCorrespondingNominalAmp,
      const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
          &ismCorrectionCalzetti) const;
  std::vector<CLineRatioCatalog> m_lineRatioCatalogs;
};

} // namespace NSEpic

#endif

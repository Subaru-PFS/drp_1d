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
#include "RedshiftLibrary/linemodel/element.h"
#include "RedshiftLibrary/operator/pdfz.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"

#include <boost/format.hpp>

#include <string>
#include <vector>

namespace NSEpic {

/**
 * \ingroup Redshift
 */
class CLineCatalogsTplRatio {

public:
  Float64 GetBestFit(const TInt32List &validLinesIndex,
                     const TFloat64List &fittedAmplitudes,
                     const TFloat64List &fittedErrors,
                     TFloat64List &amplitudesCorrected,
                     std::string &bestTplName) const;
  CLineMap const &GetRestLinesList(Int32 index) const {
    return m_lineRatioCatalogs[index].GetList();
  };
  Int32 GetCatalogsCount() const { return m_lineRatioCatalogs.size(); }

  TFloat64List getCatalogsPriors() const;
  std::string GetCatalogName(Int32 idx) const {
    return m_lineRatioCatalogs.at(idx).getName();
  }
  Int32 GetIsmIndex(Int32 idx) const {
    return m_lineRatioCatalogs.at(idx).getIsmIndex();
  }

  bool GetCatalogVelocities(Int32 idx, Float64 &elv, Float64 &alv) const;
  std::vector<std::vector<TFloat64List>> InitLineCorrespondingAmplitudes(
      const std::vector<TLineModelElementParam_ptr>
          &LineModelElementList, // TODO refactor this, this should not be a
                                 // linemodelelementlist, because all the
                                 // information retrieved from it are not
                                 // observation/spectrum dependant
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
  Float64 computeFitValue(const TFloat64List &ampsLM, const TFloat64List &errLM,
                          const TFloat64List &ampsTPL,
                          TFloat64List &ampsCorrected) const;
  Float64 getFitForOneCatalog(const TInt32List &validLinesIndex,
                              const TFloat64List &fittedAmplitudes,
                              const TFloat64List &fittedErrors,
                              const CLineRatioCatalog &catalog,
                              TFloat64List &ampsCorrected) const;
  void logLineNominalAmp(
      const std::vector<TLineModelElementParam_ptr> &LineModelElementList,
      bool enableISMCalzetti,
      const std::vector<std::vector<TFloat64List>>
          &lineCatalogLinesCorrespondingNominalAmp,
      const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
          &ismCorrectionCalzetti) const;
  std::vector<CLineRatioCatalog> m_lineRatioCatalogs;
};

} // namespace NSEpic

#endif

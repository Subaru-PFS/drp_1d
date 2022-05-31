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
#ifndef _REDSHIFT_LINE_CATALOGSTPLSHAPE_
#define _REDSHIFT_LINE_CATALOGSTPLSHAPE_

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
class CLineCatalogsTplShape {

public:
  bool Init(Int32 enableISMCalzetti,
            const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
                &ismCorrectionCalzetti,
            Float64 nsigmasupport);

  // bool AreCatalogsAligned( const CLineCatalog::TLineVector& restLineList,
  // Int32 typeFilter, Int32 forceFilter  );
  Float64 GetBestFit(const CLineCatalog::TLineVector &restLineList,
                     const TFloat64List &fittedAmplitudes,
                     const TFloat64List &fittedErrors,
                     TFloat64List &amplitudesCorrected,
                     std::string &bestTplName) const;
  CLineCatalog::TLineVector GetRestLinesList(Int32 index) const;
  Int32 GetCatalogsCount() const;
  const TFloat64List &getCatalogsPriors();
  std::string GetCatalogName(Int32 idx) const;
  Int32 GetIsmIndex(Int32 idx) const;
  Float64 GetIsmCoeff(Int32 idx) const;

  bool GetCatalogVelocities(Int32 idx, Float64 &elv, Float64 &alv) const;
  bool InitLineCorrespondingAmplitudes(
      const CLineModelElementList &LineModelElementList);
  const CLineCatalog &GetCatalog(Int32 icatlog) const;
  const std::vector<std::vector<TFloat64List>> &
  getNominalAmplitudeCorrespondance() const {
    return m_LineCatalogLinesCorrespondingNominalAmp;
  };

  void addLineRatioCatalog(const CLineRatioCatalog &lr_catalog) {
    m_lineRatioCatalogs.push_back(lr_catalog);
  }

private:
  Float64 GetFit(const TFloat64List &ampsLM, const TFloat64List &errLM,
                 const TFloat64List &ampsTPL,
                 TFloat64List &ampsCorrected) const;

  std::vector<CLineRatioCatalog> m_lineRatioCatalogs;
  std::vector<std::vector<TFloat64List>>
      m_LineCatalogLinesCorrespondingNominalAmp;
  TFloat64List m_catalogsPriors; // TODO temporary hack before reviewing
                                 // getCatalogsPrior uses

  std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
      m_ismCorrectionCalzetti;
  Int32 m_opt_dust_calzetti;
  Float64 m_nsigmasupport;
};

} // namespace NSEpic

#endif
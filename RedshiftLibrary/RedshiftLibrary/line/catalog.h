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
#ifndef _REDSHIFT_LINE_CATALOG_
#define _REDSHIFT_LINE_CATALOG_

#include <string>
#include <vector>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/line/linedetected.h"

namespace NSEpic {
class CSpectrumFluxCorrectionMeiksin;
/**
 * \ingroup Redshift
 * Line catalog allow to store multiple lines description in a single text file.
 *
 */

template <typename TLine = CLine> class CLineCatalogBase {

  static_assert(std::is_base_of<CLine, TLine>::value);

public:
  using TLineVector = std::vector<TLine>;
  using TLineMap = std::map<Int32, TLine>;

  CLineCatalogBase(Float64 sigmaSupport) : m_nSigmaSupport(sigmaSupport){};
  virtual ~CLineCatalogBase() = default;
  CLineCatalogBase() = default;
  CLineCatalogBase(const CLineCatalogBase &other) = default;
  CLineCatalogBase(CLineCatalogBase &&other) = default;
  CLineCatalogBase &operator=(const CLineCatalogBase &other) = default;
  CLineCatalogBase &operator=(CLineCatalogBase &&other) = default;

  void Add(const TLine &r) { m_List[r.GetID()] = r; };
  void Add(TLine &&r) { m_List[r.GetID()] = std::move(r); };

  const TLineMap &GetList() const { return m_List; };
  TLineMap
  GetFilteredList(CLine::EType typeFilter = CLine::EType::nType_All,
                  CLine::EForce forceFilter = CLine::EForce::nForce_All) const;
  TLineMap GetFilteredList(const std::string &typeFilter,
                           const std::string &forceFilter) const;
  static std::map<std::string, TLineVector>
  ConvertToGroupList(const TLineMap &filteredList);

  void setAsymProfileAndParams(const std::string &profile, TAsymParams params);

  void setLineAmplitude(Int32 id, Float64 nominalAmplitude);
  void convertLineProfiles2SYMIGM(
      const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrection);

protected:
  TLineMap m_List;

  Float64 m_nSigmaSupport = N_SIGMA_SUPPORT;
};

// specialization: AddLineFromParams is only implemented for
// CLineCatalogBase<CLine>
class CLineCatalog : public CLineCatalogBase<> {
public:
  using CLineCatalogBase<>::CLineCatalogBase;
  virtual ~CLineCatalog() = default;
  CLineCatalog(const CLineCatalog &other) = default;
  CLineCatalog(CLineCatalog &&other) = default;
  CLineCatalog &operator=(const CLineCatalog &other) = default;
  CLineCatalog &operator=(CLineCatalog &&other) = default;

  void AddLineFromParams(
      const std::string &name, Float64 position, std::string const &type,
      std::string const &force, const std::string &profile,
      const TAsymParams &asymParams, const std::string &groupName,
      Float64 nominalAmplitude, const std::string &velocityGroup,
      Float64 velocityOffset, bool enableVelocityFit, Int32 id,
      const std::string &str_id,
      const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrection =
          nullptr);
};

class CLineDetectedCatalog : public CLineCatalogBase<CLineDetected> {
public:
  using CLineCatalogBase<CLineDetected>::CLineCatalogBase;
  void Add(const CLineDetected &r, Int32 id) = delete;
  void Add(const CLineDetected &r) { m_List[m_List.size()] = r; };

  void Sort();
};

template class CLineCatalogBase<>;
template class CLineCatalogBase<CLineDetected>;

} // namespace NSEpic

#endif

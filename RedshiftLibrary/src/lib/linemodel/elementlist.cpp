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
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/linemodelsolution.h"

#include <cfloat>
#include <cmath>
#include <numeric>

using namespace NSEpic;
using namespace std;
/**
 * \brief Returns the number of m_Elements that fail IsOutsideLambdaRange().
 **/
Int32 CLineModelElementList::GetModelValidElementsNDdl() const {
  Int32 nddl = 0;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true) {
      continue;
    }

    nddl++;
  }
  return nddl;
}

/**
 * \brief Returns the number of elements that have only subelements with
 *non-positive amplitude.
 **/
Int32 CLineModelElementList::GetModelNonZeroElementsNDdl() const {
  Int32 nddl = 0;
  for (auto const &elt : m_Elements) {
    if (elt->IsOutsideLambdaRange())
      continue;
    if (!elt->isAllAmplitudesNull())
      nddl++;
  }
  return nddl;
}

/**
 * \brief Returns the list of indexes of elements that fail
 *IsOutsideLambdaRange.
 **/
TInt32List CLineModelElementList::GetModelValidElementsIndexes() const {
  TInt32List nonZeroIndexes;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true)
      continue;

    if (IsElementIndexInDisabledList(iElts))
      continue;

    nonZeroIndexes.push_back(iElts);
  }
  return nonZeroIndexes;
}

bool CLineModelElementList::IsElementIndexInDisabledList(Int32 index) const {
  return std::find(m_elementsDisabledIndexes.cbegin(),
                   m_elementsDisabledIndexes.cend(),
                   index) != m_elementsDisabledIndexes.cend();
}

/**
 * @brief CLineModelElementList::SetElementIndexesDisabledAuto
 * Disables all the elements that have all sub-elements (lines) amplitudes equal
 * to zero
 */
void CLineModelElementList::SetElementIndexesDisabledAuto() {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    auto const &elt = m_Elements[iElts];
    if (elt->IsOutsideLambdaRange())
      continue;

    if (elt->isAllAmplitudesNull())
      m_elementsDisabledIndexes.push_back(iElts);
  }
}

void CLineModelElementList::ResetElementIndexesDisabled() {
  m_elementsDisabledIndexes.clear();
}

/**
 * \brief Returns the list of groups, with each group being a set of line
 *indexes with the velcocity to be jointly TEMPORARY-DEV: return all the
 *indexes individually as  agroup
 **/
std::vector<TInt32List>
CLineModelElementList::GetModelVelfitGroups(CLine::EType lineType) const {
  TStringList tags;
  TInt32List nonGroupedLines;

  TInt32List nonZeroIndexes = GetModelValidElementsIndexes();

  std::map<std::string, TInt32List> tag_groups;
  for (auto iElts : nonZeroIndexes) {
    for (const auto &line : m_Elements[iElts]->GetLines())
      if (lineType == line.GetType()) {
        std::string tag = line.GetVelGroupName();
        if (tag == undefStr)
          tag = "single_" + line.GetStrID();
        tag_groups[tag].push_back(iElts);
      }
  }

  // remove duplicate elements in each groups
  for (auto [_, group] : tag_groups) {
    std::sort(group.begin(), group.end());
    group.erase(std::unique(group.begin(), group.end()), group.end());
  }

  // print the groups
  for (auto const &[tag, group] : tag_groups) {
    Log.LogDebug("    model: Group %s: nlines=%d", tag.c_str(), group.size());
    for (auto const &element : group) {
      Log.LogDebug("    model: \t iElt=%d", element);
    }
  }

  // return vector, removing tag keys
  std::vector<TInt32List> groups;
  groups.reserve(tag_groups.size());
  for (auto [_, group] : tag_groups)
    groups.push_back(std::move(group));

  return groups;
}

/**
 * \brief Returns a sorted, de-duplicated list of indices of lines whose support
 *overlap ind's support and are not listed in the argument excludedInd.
 **/

TInt32List CLineModelElementList::getOverlappingElements(
    Int32 ind, const TInt32Set &excludedInd, Float64 redshift,
    Float64 overlapThres) const {
  TInt32List indexes;

  const auto &refElement = *m_Elements[ind];

  if (refElement.IsOutsideLambdaRange()) {
    indexes.push_back(ind);
    return indexes;
  }

  auto const &refLinesList = refElement.GetLines();
  auto const refLineType = refElement.GetElementType();

  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    const auto &element = *m_Elements[iElts];

    // skip itself
    if (iElts == ind) {
      indexes.push_back(ind);
      continue;
    }

    if (element.GetElementType() != refLineType)
      continue;

    if (element.IsOutsideLambdaRange())
      continue;

    // check if in exclusion list
    if (excludedInd.find(iElts) != excludedInd.cend())
      continue;

    auto const &linesElt = element.GetLines();

    for (Int32 eltLineIdx = 0; eltLineIdx != element.GetSize(); ++eltLineIdx) {
      auto eltLine = element.GetLines()[eltLineIdx];
      if (element.IsOutsideLambdaRange(eltLineIdx))
        continue;
      for (Int32 refLineIdx = 0; refLineIdx != refLinesList.size();
           ++refLineIdx) {
        auto const &refLine = refLinesList[refLineIdx];
        if (refElement.IsOutsideLambdaRange(refLineIdx))
          continue;
        const Float64 muRef = refLine.GetPosition() * (1 + redshift);
        const Float64 sigmaRef =
            refElement.GetLineWidth(muRef, refLine.IsEmission());
        const Float64 winsizeRef =
            refLine.GetProfile()->GetNSigmaSupport() * sigmaRef;
        const Float64 overlapSizeMin = winsizeRef * overlapThres;
        const Float64 xinf = muRef - winsizeRef / 2.0;
        const Float64 xsup = muRef + winsizeRef / 2.0;

        const Float64 muElt = eltLine.GetPosition() * (1 + redshift);
        const Float64 sigmaElt =
            element.GetLineWidth(muElt, eltLine.IsEmission());
        const Float64 winsizeElt =
            eltLine.GetProfile()->GetNSigmaSupport() * sigmaElt;
        const Float64 yinf = muElt - winsizeElt / 2.0;
        const Float64 ysup = muElt + winsizeElt / 2.0;

        const Float64 max = std::max(xinf, yinf);
        const Float64 min = std::min(xsup, ysup);
        if (max - min < -overlapSizeMin) {
          indexes.push_back(iElts);
          break;
        }
      }
    }
  }

  std::sort(indexes.begin(), indexes.end());
  indexes.erase(std::unique(indexes.begin(), indexes.end()), indexes.end());

  return indexes;
}

/**
 * \brief If argument j is a valid index of m_Elements, updates the element in
 *that index calling its SetFittedAmplitude with arguments a and snr.
 **/
void CLineModelElementList::SetElementAmplitude(Int32 j, Float64 a,
                                                Float64 snr) {
  if (j >= 0 && j < m_Elements.size()) {
    m_Elements[j]->SetElementAmplitude(a, snr);
  }
  return;
}

/**
 * \brief If j is a valid index of m_Elements, returns a call to that element's
 *GetElementAmplitude. If not, returns -1.
 **/
Float64 CLineModelElementList::GetElementAmplitude(Int32 j) const {
  Float64 a = -1.0;
  if (j >= 0 && j < m_Elements.size()) {
    a = m_Elements[j]->GetElementAmplitude();
  }
  return a;
}

/**
 * @brief: get all valid elements based on whether the amplitude is null
 *
 * @param lineTypeFilter
 * @return TFloat64List
 */
TInt32List CLineModelElementList::getValidElementIndices(
    CLine::EType lineTypeFilter) const {

  const TInt32List validEltsIdx = GetModelValidElementsIndexes();
  TInt32List nonZeroValidEltsIdx;
  for (const Int32 eIdx : validEltsIdx) {
    auto const lineType = m_Elements[eIdx]->GetElementType();
    if (lineTypeFilter != CLine::EType::nType_All && lineTypeFilter != lineType)
      continue;

    if (!std::isnan(m_Elements[eIdx]->GetElementAmplitude()) &&
        m_Elements[eIdx]->GetElementAmplitude() > 0.0)
      nonZeroValidEltsIdx.push_back(eIdx);
  }
  return nonZeroValidEltsIdx;
}

TInt32List
CLineModelElementList::findElementTypeIndices(CLine::EType type) const {
  TInt32List idx;
  for (size_t i = 0; i < m_Elements.size(); ++i)
    if (m_Elements[i]->GetElementType() == type)
      idx.push_back(i);
  return idx;
}

/**
 * \brief Returns the first index of m_Elements where calling the element's
 *findElementIndex method with LineCatalogIndex argument does not return -1.
 * Returns also the line index
 **/
std::pair<Int32, Int32>
CLineModelElementList::findElementIndex(Int32 line_id) const {
  Int32 elt_index = undefIdx;
  Int32 line_index = undefIdx;
  auto it = std::find_if(begin(), end(), [line_id](auto const &elt_ptr) {
    return elt_ptr->hasLine(line_id);
  });
  if (it != end()) {
    elt_index = it - begin();
    line_index = it->get()->getLineIndex(line_id);
  }
  return std::make_pair(elt_index, line_index);
}

std::pair<Int32, Int32>
CLineModelElementList::findElementIndex(const std::string &LineTagStr,
                                        CLine::EType linetype) const {
  Int32 elt_index = undefIdx;
  Int32 line_index = undefIdx;
  auto it = std::find_if(begin(), end(),
                         [&LineTagStr, &linetype](auto const &elt_ptr) {
                           return elt_ptr->hasLine(LineTagStr) &&
                                  (linetype == CLine::EType::nType_All ||
                                   elt_ptr->GetElementType() == linetype);
                         });
  if (it != end()) {
    elt_index = it - begin();
    line_index = it->get()->getLineIndex(LineTagStr);
  }
  return std::make_pair(elt_index, line_index);
}

// first element returned is the list of element indices, then all remaining
// elts are the list of line indices, for each element listed.
std::vector<std::pair<Int32, TInt32List>>
CLineModelElementList::getIgmLinesIndices() const {

  std::vector<std::pair<Int32, TInt32List>> indices;

  for (Int32 elt_idx = 0; elt_idx < m_Elements.size(); ++elt_idx) {
    TInt32List const &line_indices = m_Elements[elt_idx]->getIgmLinesIndices();
    if (!line_indices.empty())
      indices.push_back({elt_idx, line_indices});
  }
  return indices;
}

/**
 * \brief Returns a sorted set of samples indices present in the supports of the
 *argument. For each EltsIdx entry, if the entry is not outside lambda range,
 *get the support of each subelement. For each selected support, get the sample
 *index. Sort this list and remove multiple entries. Return this clean list.
 **/
TInt32List
CLineModelElementList::getSupportIndexes(const TInt32List &EltsIdx) const {
  TInt32List indexes;

  if (!EltsIdx.size())
    return indexes;
  TInt32RangeList support;
  for (Int32 iElts : EltsIdx) {
    if (m_Elements[iElts]->IsOutsideLambdaRange())
      continue;

    TInt32RangeList s = m_Elements[iElts]->getSupport();
    support.insert(support.end(), s.begin(), s.end());
  }

  for (Int32 iS = 0; iS < support.size(); iS++) {
    for (Int32 j = support[iS].GetBegin(); j <= support[iS].GetEnd(); j++)
      indexes.push_back(j);
  }

  std::sort(indexes.begin(), indexes.end());
  indexes.erase(std::unique(indexes.begin(), indexes.end()), indexes.end());

  return indexes;
}

void CLineModelElementList::addToSpectrumAmplitudeOffset(
    const CSpectrumSpectralAxis &spectralAxis, CSpectrumFluxAxis &modelfluxAxis,
    const TInt32List &eIdx_list, CLine::EType lineTypeFilter) const {

  const auto ampOffsetGroups = getFittingGroups(eIdx_list, lineTypeFilter);

  Log.LogDetail("Elementlist: Adding n=%d ampOffsets", ampOffsetGroups.size());

  // need to avoid overlapping polynomes since the ovelapping condition is
  // computed unsing the lambda range support of the lines that intersect with a
  // bigger fraction than OVERLAP_THRES_HYBRID_FIT. It can lead to shared pixels
  // between elements considered not overlapped. In this case we should not add
  // the overlapped polynomes, but choose one of them.
  TInt32List mask(modelfluxAxis.GetSamplesCount(), 1);
  for (const auto &g : ampOffsetGroups) {
    const auto &group_eIdx_list = g.second;
    auto samples = getSupportIndexes(group_eIdx_list);
    const auto &pCoeffs =
        m_Elements[group_eIdx_list.front()]->GetPolynomCoeffs();
    for (Int32 s : samples) {
      modelfluxAxis[s] += pCoeffs.getValue(spectralAxis[s]) * mask[s];
      mask[s] = 0; // one sample can only be written once
    }
  }
}

void CLineModelElementList::resetAmplitudeOffset() {
  for (auto &elt : m_Elements) {
    elt->SetPolynomCoeffs(TPolynomCoeffs());
  }
}

/**
 * \brief Get the scale marginalization correction
 *
 * WARNING: (todo-check) for lm-rules or lm-free, if hybrid method was used to
 *fit, mtm and dtm are not estimated for now...
 **/
Float64 CLineModelElementList::getScaleMargCorrection(Int32 Eltidx) const {
  Float64 corr = 0.0;

  // scale marg for continuum
  // corr += getContinuumScaleMargCorrection();
  Int32 iElts_start = 0;
  Int32 iElts_end = m_Elements.size();
  if (Eltidx != undefIdx) {
    iElts_start = Eltidx;
    iElts_end = Eltidx + 1;
  }
  for (Int32 iElts = iElts_start; iElts != iElts_end; iElts++) {
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true)
      continue;

    Float64 mtm = m_Elements[iElts]->GetSumGauss();
    if (mtm > 0.0)
      corr += log(mtm);
  }

  return corr;
}

/**
 * @brief CLineModelFitting::GetModelStrongLinePresent
 * @return 1 if there is 1 strong emission line present
 */
bool CLineModelElementList::GetModelStrongEmissionLinePresent() const {

  TInt32List validEltsIdx = GetModelValidElementsIndexes();
  for (Int32 iElts : validEltsIdx) {
    auto const &elt = m_Elements[iElts];
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      auto const &line = elt->GetLines()[index];
      if (!line.IsEmission() || !line.IsStrong())
        continue;

      Float64 amp = elt->GetFittedAmplitude(index);
      if (amp > 0.0) {
        Log.LogDebug("    model: GetModelStrongEmissionLinePresent - found "
                     "Strong EL: %s",
                     line.GetName().c_str());
        return true;
      }
    }
  }

  return false;
}

/**
 * @brief
 * @return 1 if ha em is the strongest line present
 */
bool CLineModelElementList::GetModelHaStrongest() const {

  Float64 ampMax = -DBL_MAX;
  std::string ampMaxLineTag = "";

  TInt32List validEltsIdx = GetModelValidElementsIndexes();
  for (Int32 iElts : validEltsIdx) {
    auto const &elt = m_Elements[iElts];
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      auto const &line = elt->GetLines()[index];
      if (!line.IsEmission())
        continue;

      Float64 amp = elt->GetFittedAmplitude(index);
      if (amp > 0. && amp > ampMax) {
        ampMaxLineTag = line.GetName().c_str();
        ampMax = amp;
      }
    }
  }

  bool isHaStrongest = (!std::isnan(ampMax) && ampMax > 0. &&
                        ampMaxLineTag == linetags::halpha_em);
  if (isHaStrongest) {
    Log.LogDebug("    model: GetModelHaStrongest - found to be true with "
                 "ampMax=%e (for line=Halpha)",
                 ampMax);
  }
  return isHaStrongest;
}

std::map<std::string, TInt32List>
CLineModelElementList::getFittingGroups(TInt32List EltsIdx,
                                        CLine::EType lineTypeFilter) const {
  std::map<std::string, TInt32List> fittingGroups;

  // if not elments list do with all
  if (EltsIdx.empty()) {
    EltsIdx.resize(size());
    std::iota(EltsIdx.begin(), EltsIdx.end(), 0);
  }

  for (auto idx : EltsIdx) {
    const auto &elt = m_Elements[idx];

    auto const lineType = elt->GetElementType();
    if (lineTypeFilter != CLine::EType::nType_All && lineTypeFilter != lineType)
      continue;

    const auto &info = elt->GetFittingGroupInfo();
    if (info == undefStr)
      continue;
    if (fittingGroups.find(info) == fittingGroups.end()) {
      fittingGroups[info] = {idx};
      continue;
    }

    fittingGroups[info].push_back(idx);
  }

  return fittingGroups;
}

/**
 * @brief Look for polynom coeffs corresponding to one specific Line
 *
 * @param eIdx
 * @return TPolynomCoeffs
 */
const TPolynomCoeffs &
CLineModelElementList::getPolynomCoeffs(Int32 eIdx) const {
  return m_Elements[eIdx]->GetPolynomCoeffs();
}

/**
 * @brief
 *
 * @param eIdx_list
 * @param subeIdx_list
 * @param sigma_support
 * @return TInt32RangeList
 */
TInt32RangeList CLineModelElementList::getlambdaIndexesUnderLines(
    const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
    Float64 sigma_support, const CSpectrumSpectralAxis &spectralAxis,
    const TFloat64Range &lambdaRange, Float64 redshift) const {

  TInt32RangeList indexRangeList(eIdx_list.size());
  for (Int32 i = 0; i < eIdx_list.size(); i++) {
    Int32 eIdx = eIdx_list[i];
    Int32 subeIdx = subeIdx_list[i];

    Float64 mu = NAN;
    Float64 LineWidth = NAN;
    m_Elements[eIdx]->getObservedPositionAndLineWidth(subeIdx, redshift, mu,
                                                      LineWidth);

    Float64 winsizeAngstrom = LineWidth * sigma_support;

    indexRangeList[i] = CLineModelElement::EstimateIndexRange(
        spectralAxis, mu, lambdaRange, winsizeAngstrom);
  }

  if (eIdx_list.size() == 1)
    return indexRangeList;
  TInt32RangeList nonOverlappingIndexRangeList =
      TInt32Range::joinIntersections(std::move(indexRangeList));

  return nonOverlappingIndexRangeList;
}

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
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true) {
      continue;
    }
    bool isAllZero = true;
    for (Int32 ie = 0; ie < m_Elements[iElts]->GetSize(); ie++) {
      if (m_Elements[iElts]->GetFittedAmplitude(ie) > 0.0) {
        isAllZero = false;
      }
    }

    if (isAllZero == false) {
      nddl++;
    }
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
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true) {
      continue;
    }
    if (IsElementIndexInDisabledList(iElts)) {
      continue;
    }

    nonZeroIndexes.push_back(iElts);
  }
  return nonZeroIndexes;
}

bool CLineModelElementList::IsElementIndexInDisabledList(Int32 index) const {
  for (Int32 i = 0; i < m_elementsDisabledIndexes.size(); i++) {
    if (m_elementsDisabledIndexes[i] == index) {
      return true;
    }
  }
  return false;
}

/**
 * @brief CLineModelElementList::SetElementIndexesDisabledAuto
 * Disables all the elements that have all sub-elements (lines) amplitudes equal
 * to zero
 */
void CLineModelElementList::SetElementIndexesDisabledAuto() {
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true) {
      continue;
    }
    bool isAllZero = true;
    for (Int32 ie = 0; ie < m_Elements[iElts]->GetSize(); ie++) {
      if (m_Elements[iElts]->GetFittedAmplitude(ie) > 0.0) {
        isAllZero = false;
      }
    }

    if (isAllZero == true) {
      m_elementsDisabledIndexes.push_back(iElts);
    }
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
CLineModelElementList::GetModelVelfitGroups(Int32 lineType) const {
  TStringList tags;
  TInt32List nonGroupedLines;

  TInt32List nonZeroIndexes = GetModelValidElementsIndexes();
  for (auto iElts : nonZeroIndexes)
    for (const auto &line : m_Elements[iElts]->GetLines())
      if (lineType == line.GetType()) {
        std::string _tag = line.GetVelGroupName();
        if (_tag != "-1")
          tags.push_back(_tag);
        else
          nonGroupedLines.push_back(iElts);
      }

  // create the group tag set by removing duplicates
  std::sort(tags.begin(), tags.end());
  tags.erase(std::unique(tags.begin(), tags.end()), tags.end());

  // add the grouped lines
  std::vector<TInt32List> groups;
  TStringList groupsTags;
  for (Int32 itag = 0; itag < tags.size(); itag++) {
    TInt32List _group;
    for (auto iElts : nonZeroIndexes)
      for (const auto &line : m_Elements[iElts]->GetLines())
        if (lineType == line.GetType())
          if (tags[itag] == line.GetVelGroupName())
            _group.push_back(iElts);

    // add the grouped lines, no duplicates
    std::sort(_group.begin(), _group.end());
    _group.erase(std::unique(_group.begin(), _group.end()), _group.end());

    groups.push_back(_group);
    groupsTags.push_back(tags[itag]);
  }
  // add the non grouped lines, no duplicates
  std::sort(nonGroupedLines.begin(), nonGroupedLines.end());
  nonGroupedLines.erase(
      std::unique(nonGroupedLines.begin(), nonGroupedLines.end()),
      nonGroupedLines.end());
  for (const auto lIdx : nonGroupedLines) {
    groups.push_back({lIdx});
    groupsTags.push_back("-1");
  }

  if (true) {
    // print the groups
    for (Int32 igr = 0; igr < groups.size(); igr++) {
      Log.LogDebug("    model: Group %d/%d: nlines=%d, tag=%s", igr + 1,
                   groups.size(), groups[igr].size(), groupsTags[igr].c_str());
      for (Int32 i = 0; i < groups[igr].size(); i++) {
        Log.LogDebug("    model: \t%d: iElt=%d", i + 1, groups[igr][i]);
      }
    }
  }

  // Override velGroups from Catalog: => Individual lines as groups
  /*
  std::vector<TInt32List> groups;
  TInt32List nonZeroIndexes = GetModelValidElementsIndexes();
  for(Int32 i=0; i<nonZeroIndexes.size(); i++)
  {
      if(lineType == m_Elements[nonZeroIndexes[i]]->GetElementType())
      {
          TInt32List gr;
          gr.push_back(nonZeroIndexes[i]);
          groups.push_back(gr);
          //Log.LogInfo("Group %d, idx=%d", groups.size(),
  groups[groups.size()-1][0]);
      }
  }
  //*/

  return groups;
}

/**
 * \brief Returns a sorted, de-duplicated list of indices of lines whose support
 *overlap ind's support and are not listed in the argument excludedInd.
 **/

TInt32List CLineModelElementList::getOverlappingElements(
    Int32 ind, const TInt32List &excludedInd, Float64 redshift,
    Float64 overlapThres) const {
  TInt32List indexes;

  if (m_Elements[ind]->IsOutsideLambdaRange()) {
    indexes.push_back(ind);
    return indexes;
  }

  const TLineVector &linesRef = m_Elements[ind]->GetLines();
  Int32 linetypeRef = m_Elements[ind]->GetElementType();

  Int32 xinf = 0;
  Int32 yinf = 0;
  Int32 xsup = 0;
  Int32 ysup = 0;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    // check linetype
    if (m_Elements[iElts]->GetElementType() != linetypeRef) {
      continue;
    }

    // check if outside lambdarange
    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      continue;
    }

    // check if in exclusion list
    bool excluded = false;
    for (Int32 iexcl = 0; iexcl < excludedInd.size(); iexcl++) {
      if (iElts == excludedInd[iexcl]) {
        excluded = true;
        break;
      }
    }
    if (excluded) {
      continue;
    }

    const TLineVector &linesElt = m_Elements[iElts]->GetLines();

    for (Int32 iLineElt = 0; iLineElt < linesElt.size(); iLineElt++) {
      for (Int32 iLineRef = 0; iLineRef < linesRef.size(); iLineRef++) {
        Float64 muRef = linesRef[iLineRef].GetPosition() * (1 + redshift);
        Float64 cRef = m_Elements[ind]->GetLineWidth(
            muRef, linesRef[iLineRef].GetIsEmission());
        Float64 winsizeRef =
            linesRef[iLineRef].GetProfile()->GetNSigmaSupport() * cRef;
        Float64 overlapSizeMin = winsizeRef * overlapThres;
        xinf = muRef - winsizeRef / 2.0;
        xsup = muRef + winsizeRef / 2.0;

        Float64 muElt = linesElt[iLineElt].GetPosition() * (1 + redshift);
        Float64 cElt = m_Elements[iElts]->GetLineWidth(
            muElt, linesElt[iLineElt].GetIsEmission());
        Float64 winsizeElt =
            linesElt[iLineElt].GetProfile()->GetNSigmaSupport() * cElt;
        yinf = muElt - winsizeElt / 2.0;
        ysup = muElt + winsizeElt / 2.0;

        Float64 max = std::max(xinf, yinf);
        Float64 min = std::min(xsup, ysup);
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

// TODO never called function, keep it ?
/**
 * \brief Returns a sorted set of line indices present in the supports of the
 *argument. Create a vector named indexes. If the argument ind is an index to
 *m_Elements that IsOutSideLambdaRange, return indexes. For each entry in
 *m_Elements: If the entry has a different linetype than the line
 *corresponding to ind, go to the next entry. If the entry
 *IsOutsideLambdaRange, go to the enxt entry. For each subentry in the support
 *of entry: For each subsubentry in the support of ind: If the overlap in the
 *spectralAxis is smaller than -1 * overlapThres * winsize, add entry to
 *indexes. Sort indexes, remove duplicates from indexes, and return indexes.
 **/
TInt32List CLineModelElementList::getOverlappingElementsBySupport(
    Int32 ind, Float64 redshift, const CSpectrumSpectralAxis &spectralAxis,
    Float64 overlapThres) const {

  TInt32List indexes;

  if (m_Elements[ind]->IsOutsideLambdaRange()) {
    indexes.push_back(ind);
    return indexes;
  }
  TInt32RangeList refsupport = m_Elements[ind]->getSupport();
  const CLine &line = m_Elements[ind]->GetLines().front();
  Int32 linetype = line.GetType();
  Float64 mu = line.GetPosition() * (1 + redshift);
  Float64 c = m_Elements[ind]->GetLineWidth(mu, line.GetIsEmission());
  Float64 winsize = line.GetProfile()->GetNSigmaSupport() * c;
  Float64 overlapThresholdMin = winsize * overlapThres;
  // overlapThresholdMin = 0.0;

  Int32 x1 = 0;
  Int32 y1 = 0;
  Int32 x2 = 0;
  Int32 y2 = 0;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (m_Elements[iElts]->GetElementType() != linetype) {
      continue;
    }

    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      continue;
    }
    TInt32RangeList s = m_Elements[iElts]->getSupport();
    for (Int32 iS = 0; iS < s.size(); iS++) {
      for (Int32 iRefS = 0; iRefS < refsupport.size(); iRefS++) {
        x1 = refsupport[iRefS].GetBegin();
        x2 = refsupport[iRefS].GetEnd();
        y1 = s[iS].GetBegin();
        y2 = s[iS].GetEnd();

        // Log.LogInfo( "hybrid fit: iRefS=%d - support=%d,%d", iRefS, x1,
        // x2); Log.LogInfo( "hybrid fit: iS=%d - support=%d,%d", iS, y1, y2);

        //                if( std::max(x1,y1) < std::min(x2,y2) ){
        //                    indexes.push_back(iElts);
        //                    break;
        //                }

        Float64 max = spectralAxis[std::max(x1, y1)];
        Float64 min = spectralAxis[std::min(x2, y2)];
        if (max - min < -overlapThresholdMin) {
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
TInt32List
CLineModelElementList::getValidElementIndices(Int32 lineTypeFilter) const {

  const TInt32List validEltsIdx = GetModelValidElementsIndexes();
  TInt32List nonZeroValidEltsIdx;
  for (const Int32 eIdx : validEltsIdx) {
    const Int32 lineType = m_Elements[eIdx]->GetElementType();
    if (lineTypeFilter != -1 && lineTypeFilter != lineType)
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
Int32 CLineModelElementList::findElementIndex(Int32 LineCatalogIndex,
                                              Int32 &lineIdx) const {
  Int32 idx = undefIdx;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    lineIdx = m_Elements[iElts]->findElementIndex(LineCatalogIndex);
    if (lineIdx != undefIdx) {
      idx = iElts;
      break;
    }
  }
  return idx;
}

/**
 * \brief Returns the first index of m_Elements where calling the element's
 *findElementIndex method with LineTagStr argument does not return -1.
 **/
Int32 CLineModelElementList::findElementIndex(const std::string &LineTagStr,
                                              Int32 linetype,
                                              Int32 &lineIdx) const {
  Int32 idx = undefIdx;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    lineIdx = m_Elements[iElts]->findElementIndex(LineTagStr);
    if (lineIdx == undefIdx)
      continue;

    if (linetype != -1 && m_Elements[iElts]->GetElementType() != linetype)
      continue;

    idx = iElts;
    break;
  }
  return idx;
}

Int32 CLineModelElementList::findElementIndex(const std::string &LineTagStr,
                                              Int32 linetype) const {
  Int32 lineIdx = undefIdx;
  return findElementIndex(LineTagStr, linetype, lineIdx);
}

// first element returned is the list of element indices, then all remaining
// elts are the list of line indices, for each element listed.
std::tuple<TInt32List, std::vector<TInt32List>>
CLineModelElementList::getIgmLinesIndices() const {
  TInt32List EltsIdxList;
  std::vector<TInt32List> lineIdxList;
  for (Int32 iElts = 0; iElts < m_Elements.size(); ++iElts) {
    TInt32List lineIdx = m_Elements[iElts]->getIgmLinesIndices();
    if (!lineIdx.empty()) {
      EltsIdxList.push_back(iElts);
      lineIdxList.push_back(std::move(lineIdx));
    }
  }
  return std::make_tuple(std::move(EltsIdxList), std::move(lineIdxList));
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
  for (Int32 i = 0; i < EltsIdx.size(); i++) {
    Int32 iElts = EltsIdx[i];

    if (m_Elements[iElts]->IsOutsideLambdaRange()) {
      continue;
    }
    TInt32RangeList s = m_Elements[iElts]->getSupport();
    support.insert(support.end(), s.begin(), s.end());
  }

  for (Int32 iS = 0; iS < support.size(); iS++) {
    for (Int32 j = support[iS].GetBegin(); j <= support[iS].GetEnd(); j++) {
      indexes.push_back(j);
    }
  }

  std::sort(indexes.begin(), indexes.end());
  indexes.erase(std::unique(indexes.begin(), indexes.end()), indexes.end());

  return indexes;
}

void CLineModelElementList::addToSpectrumAmplitudeOffset(
    const CSpectrumSpectralAxis &spectralAxis, CSpectrumFluxAxis &modelfluxAxis,
    const TInt32List &eIdx_list, Int32 lineTypeFilter) const {

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
Float64 CLineModelElementList::getScaleMargCorrection(Int32 idxLine) const {
  Float64 corr = 0.0;

  // scale marg for continuum
  // corr += getContinuumScaleMargCorrection();

  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (idxLine != undefIdx && idxLine != iElts) {
      continue;
    }
    if (m_Elements[iElts]->IsOutsideLambdaRange() == true) {
      continue;
    }
    // if(m_Elements[iElts]->GetElementAmplitude()<=0.0){
    //     continue;
    // }

    Float64 mtm = m_Elements[iElts]->GetSumGauss();
    if (mtm > 0.0) {
      corr += log(mtm);
    }
  }

  return corr;
}

/**
 * @brief CLineModelFitting::GetModelStrongLinePresent
 * @return 1 if there is 1 strong emission line present
 */
bool CLineModelElementList::GetModelStrongEmissionLinePresent() const {

  // TODO is this check really necessary here ?
  //   if (!m_RestLineList.size())
  //   THROWG(INTERNAL_ERROR, "m_RestframeList is empty");

  bool isStrongPresent = false;

  TInt32List validEltsIdx = GetModelValidElementsIndexes();
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    Int32 nlines = m_Elements[iElts]->GetSize();
    for (Int32 lineIdx = 0; lineIdx < nlines; lineIdx++) {
      if (!m_Elements[iElts]->GetLines()[lineIdx].GetIsEmission() ||
          !m_Elements[iElts]->GetLines()[lineIdx].GetIsStrong()) {
        continue;
      }

      Float64 amp = m_Elements[iElts]->GetFittedAmplitude(lineIdx);
      if (amp > 0.0) {
        isStrongPresent = true;
        Log.LogDebug("    model: GetModelStrongEmissionLinePresent - found "
                     "Strong EL: %s",
                     m_Elements[iElts]->GetLines()[lineIdx].GetName().c_str());
        break;
      }
    }
    if (isStrongPresent)
      break;
  }

  return isStrongPresent;
}

/**
 * @brief
 * @return 1 if ha em is the strongest line present
 */
bool CLineModelElementList::GetModelHaStrongest() const {

  Float64 ampMax = -DBL_MAX;
  std::string ampMaxLineTag = "";

  TInt32List validEltsIdx = GetModelValidElementsIndexes();
  for (Int32 iValidElts = 0; iValidElts < validEltsIdx.size(); iValidElts++) {
    Int32 iElts = validEltsIdx[iValidElts];
    Int32 nlines = m_Elements[iElts]->GetSize();
    for (Int32 lineIdx = 0; lineIdx < nlines; lineIdx++) {
      if (!m_Elements[iElts]->GetLines()[lineIdx].GetIsEmission()) {
        continue;
      }

      Float64 amp = m_Elements[iElts]->GetFittedAmplitude(lineIdx);
      if (amp > 0. && amp > ampMax) {
        ampMaxLineTag =
            m_Elements[iElts]->GetLines()[lineIdx].GetName().c_str();
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
                                        Int32 lineTypeFilter) const {
  std::map<std::string, TInt32List> fittingGroups;

  // if not elments list do with all
  if (EltsIdx.empty()) {
    EltsIdx.resize(size());
    std::iota(EltsIdx.begin(), EltsIdx.end(), 0);
  }

  for (auto idx : EltsIdx) {
    const auto &elt = m_Elements[idx];

    Int32 lineType = elt->GetElementType();
    if (lineTypeFilter != undefIdx && lineTypeFilter != lineType)
      continue;

    const auto &info = elt->GetFittingGroupInfo();
    if (info == "-1")
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

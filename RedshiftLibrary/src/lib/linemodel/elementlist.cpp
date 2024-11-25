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
#include <cfloat>
#include <cmath>
#include <numeric>

#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/linemodel/linemodelsolution.h"

using namespace NSEpic;
using namespace std;

/**
 * \brief Returns the list of indexes of elements that fail
 *IsOutsideLambdaRange.
 **/
TInt32List CLineModelElementList::GetElementsIndicesInsideLambdaRange() const {
  TInt32List validIndices;
  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    if (!m_Elements[iElts]->IsOutsideLambdaRange())
      validIndices.push_back(iElts);
  }
  return validIndices;
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

  TInt32List validIndices = GetElementsIndicesInsideLambdaRange();

  std::map<std::string, TInt32List> tag_groups;
  for (auto iElts : validIndices) {
    const auto &elt_param_ptr = m_Elements[iElts]->getElementParam();
    for (const auto &line : elt_param_ptr->GetLines())
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
    Log.LogDebug(Formatter()
                 << "    model: Group " << tag << ": nlines=" << group.size());
    for (auto const &element : group) {
      Log.LogDebug(Formatter() << "    model: \t iElt=" << element);
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
// TODO move this to LMElementListVector ( #8798 )
// la difficulté viendra avec des LSF différentes, qui donc engendre un overlap
// différent: il faudra prendre l'overlap le plus important, ie si ça overlap
// dans un spectre alors ça overlap globalement (il faut faire un fit joint)
TInt32List CLineModelElementList::getOverlappingElements(
    Int32 ind, const TInt32Set &excludedInd, Float64 redshift,
    Float64 overlapThres) const {
  TInt32List indexes;

  const auto &refElement = *m_Elements[ind];
  const auto &refElement_param = refElement.getElementParam();

  if (refElement.IsOutsideLambdaRange()) {
    indexes.push_back(ind);
    return indexes;
  }

  auto const &refLinesList = refElement_param->GetLines();
  auto const refLineType = refElement_param->GetElementType();

  for (Int32 iElts = 0; iElts < m_Elements.size(); iElts++) {
    const auto &element = *m_Elements[iElts];
    const auto &element_param = element.getElementParam();

    // skip itself
    if (iElts == ind) {
      indexes.push_back(ind);
      continue;
    }

    if (element_param->GetElementType() != refLineType)
      continue;

    if (element.IsOutsideLambdaRange())
      continue;

    // check if in exclusion list
    if (excludedInd.find(iElts) != excludedInd.cend())
      continue;

    auto const &linesElt = element_param->GetLines();

    for (Int32 eltLineIdx = 0; eltLineIdx != element.GetSize(); ++eltLineIdx) {
      auto eltLine = linesElt[eltLineIdx];
      if (element.IsOutsideLambdaRangeLine(eltLineIdx))
        continue;
      for (Int32 refLineIdx = 0; refLineIdx != refLinesList.size();
           ++refLineIdx) {
        auto const &refLine = refLinesList[refLineIdx];
        if (refElement.IsOutsideLambdaRangeLine(refLineIdx))
          continue;
        const Float64 muRef = refLine.GetPosition() * (1 + redshift);
        const Float64 sigmaRef = refElement.GetLineWidth(muRef);
        const Float64 winsizeRef =
            refLine.GetProfile()->GetNSigmaSupport() * sigmaRef;
        const Float64 overlapSizeMin = winsizeRef * overlapThres;
        const Float64 xinf = muRef - winsizeRef / 2.0;
        const Float64 xsup = muRef + winsizeRef / 2.0;

        const Float64 muElt = eltLine.GetPosition() * (1 + redshift);
        const Float64 sigmaElt = element.GetLineWidth(muElt);
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
 * @brief: get all valid elements based on whether the amplitude is null
 *
 * @param lineTypeFilter
 * @return TFloat64List
 */
TInt32List CLineModelElementList::getNonZeroElementIndices(
    CLine::EType lineTypeFilter) const {

  const TInt32List validEltsIdx = GetElementsIndicesInsideLambdaRange();
  TInt32List nonZeroValidEltsIdx;
  for (const Int32 eIdx : validEltsIdx) {
    auto const &param = m_Elements[eIdx]->getElementParam();
    auto const lineType = param->GetElementType();
    if (lineTypeFilter != CLine::EType::nType_All && lineTypeFilter != lineType)
      continue;
    Float64 const amp = param->GetElementAmplitude();
    if (!std::isnan(amp) && amp > 0.0)
      nonZeroValidEltsIdx.push_back(eIdx);
  }
  return nonZeroValidEltsIdx;
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
    return elt_ptr->getElementParam()->hasLine(line_id);
  });
  if (it != end()) {
    elt_index = it - begin();
    line_index = it->get()->getElementParam()->getLineIndex(line_id);
  }
  return std::make_pair(elt_index, line_index);
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

  Log.LogDetail(Formatter() << "Elementlist: Adding n="
                            << ampOffsetGroups.size() << " ampOffsets");

  // need to avoid overlapping polynomes since the ovelapping condition is
  // computed unsing the lambda range support of the lines that intersect with a
  // bigger fraction than OVERLAP_THRES_HYBRID_FIT. It can lead to shared pixels
  // between elements considered not overlapped. In this case we should not add
  // the overlapped polynomes, but choose one of them.
  TInt32List mask(modelfluxAxis.GetSamplesCount(), 1);
  for (const auto &[_, group_eIdx_list] : ampOffsetGroups) {
    // filter out not fittable elements
    TInt32List valid_eIdx_list;
    std::copy_if(group_eIdx_list.begin(), group_eIdx_list.end(),
                 std::back_inserter(valid_eIdx_list), [this](Int32 idx) {
                   return m_Elements[idx]->getElementParam()->isFittable();
                 });
    if (valid_eIdx_list.empty())
      continue;
    auto samples = getSupportIndexes(valid_eIdx_list);
    const auto &pCoeffs = m_Elements[valid_eIdx_list.front()]
                              ->getElementParam()
                              ->GetPolynomCoeffs();
    for (Int32 s : samples) {
      modelfluxAxis[s] += pCoeffs.getValue(spectralAxis[s]) * mask[s];
      mask[s] = 0; // one sample can only be written once
    }
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
    if (m_Elements[iElts]->getElementParam()->isNotFittable())
      continue;

    Float64 mtm = m_Elements[iElts]->getElementParam()->getSumGauss();
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

  TInt32List validEltsIdx = GetElementsIndicesInsideLambdaRange();
  for (Int32 iElts : validEltsIdx) {
    auto const &elt = m_Elements[iElts];
    auto const &elt_param = elt->getElementParam();
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      auto const &line = elt_param->GetLines()[index];
      if (!line.IsEmission() || !line.IsStrong())
        continue;

      Float64 amp = elt->getElementParam()->GetFittedAmplitude(index);
      if (amp > 0.0) {
        Log.LogDebug(Formatter()
                     << "    model: GetModelStrongEmissionLinePresent - found. "
                        "Strong EL: "
                     << line.GetName());
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

  TInt32List validEltsIdx = GetElementsIndicesInsideLambdaRange();
  for (Int32 iElts : validEltsIdx) {
    auto const &elt = m_Elements[iElts];
    auto const &elt_param = elt->getElementParam();
    for (Int32 index = 0; index != elt->GetSize(); ++index) {
      auto const &line = elt_param->GetLines()[index];
      if (!line.IsEmission())
        continue;

      Float64 amp = elt->getElementParam()->GetFittedAmplitude(index);
      if (amp > 0. && amp > ampMax) {
        ampMaxLineTag = line.GetName().c_str();
        ampMax = amp;
      }
    }
  }

  bool isHaStrongest = (!std::isnan(ampMax) && ampMax > 0. &&
                        ampMaxLineTag == linetags::halpha_em);
  if (isHaStrongest) {
    Log.LogDebug(Formatter()
                 << "    model: GetModelHaStrongest - found to be true with "
                    "ampMax="
                 << ampMax << " (for line=Halpha)");
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
    const auto &elt_param = elt->getElementParam();

    auto const lineType = elt_param->GetElementType();
    if (lineTypeFilter != CLine::EType::nType_All && lineTypeFilter != lineType)
      continue;

    const auto &info = elt_param->getFittingGroupInfo();
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
    Int32 const eIdx = eIdx_list[i];
    Int32 const subeIdx = subeIdx_list[i];
    auto const &elt = m_Elements[eIdx];
    if (sigma_support ==
        elt->getElementParam()->getLineProfile(subeIdx)->GetNSigmaSupport()) {
      // same Nsigma support, get already computed support indices
      indexRangeList[i] = elt->getTheoreticalSupportSubElt(subeIdx);
    } else {
      // different Nsigma support, thus need to recompute support indices
      auto const &[mu, LineWidth] =
          elt->getObservedPositionAndLineWidth(redshift, subeIdx);
      Float64 const winsizeAngstrom = LineWidth * sigma_support;
      indexRangeList[i] = CLineModelElement::EstimateIndexRange(
          spectralAxis, mu, lambdaRange, winsizeAngstrom);
    }
  }

  if (eIdx_list.size() == 1)
    return indexRangeList;
  TInt32RangeList const nonOverlappingIndexRangeList =
      TInt32Range::joinIntersections(std::move(indexRangeList));

  return nonOverlappingIndexRangeList;
}

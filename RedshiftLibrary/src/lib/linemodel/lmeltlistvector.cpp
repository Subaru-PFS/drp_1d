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
#include "RedshiftLibrary/linemodel/spectrummodel.h"
#include "RedshiftLibrary/processflow/autoscope.h"
#include "RedshiftLibrary/processflow/context.h"

using namespace NSEpic;
using namespace std;

CLMEltListVector::CLMEltListVector(CTLambdaRangePtrVector lambdaranges,
                                   const CSpectraGlobalIndex &spcIndex,
                                   const CLineMap &restLineList,
                                   ElementComposition element_composition)
    : m_lambdaRanges(lambdaranges), m_spectraIndex(spcIndex),
      m_RestLineList(restLineList) {

  switch (element_composition) {
  case ElementComposition::Default:
    // load the regular catalog
    LoadCatalog();
    break;
  case ElementComposition::EmissionAbsorption:
    //"tplRatio" and "tplCorr"
    // load the tplratio catalog with only 1 element for all lines
    // LoadCatalogOneMultiline(restLineList);
    // load the tplratio catalog with 2 elements: 1 for the Em lines + 1 for
    // the Abs lines
    LoadCatalogTwoMultilinesAE();
    break;
  case ElementComposition::OneLine:
    // load each line alone in one element (linemeas)
    LoadCatalogOneLineByElement();
  }
  for (auto &spcIndex : m_spectraIndex) {

    m_ElementsVector.push_back(CLineModelElementList());
    fillElements();
  }
}

CLMEltListVector::CLMEltListVector(CLineModelElementList eltlist,
                                   const CLineMap &restLineList)
    : m_RestLineList(restLineList), m_spectraIndex(1) { // for unit test

  m_ElementsVector.push_back(eltlist);

  for (auto elt : m_ElementsVector[0])
    m_ElementsParams.push_back(elt->getElementParam());
}

Int32 CLMEltListVector::GetModelNonZeroElementsNDdl() const {
  Int32 nddl = 0;
  for (Int32 elt_index = 0; elt_index < m_ElementsParams.size(); elt_index++) {
    if (m_ElementsParams[elt_index]->isNotFittable())
      continue;
    if (!m_ElementsParams[elt_index]->isAllAmplitudesNull())
      nddl++;
  }
  return nddl;
}

bool CLMEltListVector::computeOutsideLambdaRangeLine(Int32 elt_index,
                                                     Int32 line_index) {
  for (auto &spcIndex : m_spectraIndex) {
    if (!getElementList()[elt_index]->IsOutsideLambdaRangeLine(line_index))
      return false;
  }
  return true;
}

bool CLMEltListVector::computeOutsideLambdaRange(Int32 elt_index) {
  for (auto &spcIndex : m_spectraIndex) {
    if (!getElementList()[elt_index]->IsOutsideLambdaRange())
      return false;
  }
  return true;
}

void CLMEltListVector::AddElementParam(CLineVector lines) {
  size_t nb_lines = lines.size();

  CAutoScope autoscope(Context.m_ScopeStack, "lineModel");
  auto const ps = Context.GetParameterStore();
  Float64 const velocityEmission = ps->GetScoped<Float64>("velocityEmission");
  Float64 const velocityAbsorption =
      ps->GetScoped<Float64>("velocityAbsorption");
  std::string lineWidthType = ps->GetScoped<std::string>("lineWidthType");
  m_ElementsParams.push_back(std::make_shared<TLineModelElementParam>(
      std::move(lines), velocityEmission, velocityAbsorption, lineWidthType));
}

void CLMEltListVector::fillElements() {
  for (auto &ep : m_ElementsParams)
    getElementList().push_back(std::make_shared<CLineModelElement>(ep));
}
/**
 * \brief For each line in each group of the argument, finds the associated
 *line in the catalog and saves this information to getElementList(). Converts
 *the argument restLineList to a group list. For each entry in this list: For
 *each line in this entry: Finds the index in the catalog from the line name and
 *type. Saves the line, the catalog index and the nominal amplitude for the
 *line thusly associated to this line. If at least one line was found, save
 *this result in getElementList().
 **/
void CLMEltListVector::LoadCatalog() {
  auto groupList = CLineCatalog::ConvertToGroupList(m_RestLineList);
  for (auto &[_, lines] : groupList) {
    AddElementParam(std::move(lines));
  }
}

void CLMEltListVector::LoadCatalogOneLineByElement() {
  for (auto const &[_, line] : m_RestLineList) {
    AddElementParam(CLineVector{line});
  }
}

void CLMEltListVector::LoadCatalogOneMultiline() {
  CLineVector RestLineVector;
  RestLineVector.reserve(m_RestLineList.size());
  for (auto const &[_, line] : m_RestLineList)
    RestLineVector.push_back(line);

  AddElementParam(std::move(RestLineVector));
}

void CLMEltListVector::LoadCatalogTwoMultilinesAE() {

  std::vector<CLine::EType> const types = {CLine::EType::nType_Absorption,
                                           CLine::EType::nType_Emission};

  for (auto type : types) {
    CLineVector lines;
    for (auto const &[id, line] : m_RestLineList) {
      if (line.GetType() == type)
        lines.push_back(line);
    }

    if (lines.size() > 0) {
      AddElementParam(std::move(lines));
    }
  }
}

Float64 CLMEltListVector::getScaleMargCorrection(Int32 Eltidx) const {
  Float64 smc = 0;
  for (auto &spcIndex : m_spectraIndex) {

    smc += getElementList().getScaleMargCorrection(Eltidx);
  }
  return smc;
}

TInt32List CLMEltListVector::findElementTypeIndices(CLine::EType type) const {
  TInt32List idx;
  for (size_t i = 0; i < m_ElementsParams.size(); ++i)
    if (m_ElementsParams[i]->m_type == type)
      idx.push_back(i);
  return idx;
}

/**
 * \brief Returns the first index of m_Elements where calling the element's
 *findElementIndex method with LineCatalogIndex argument does not return -1.
 * Returns also the line index
 **/
std::pair<Int32, Int32>
CLMEltListVector::findElementIndex(Int32 line_id) const {
  Int32 elt_index = undefIdx;
  Int32 line_index = undefIdx;
  auto it = std::find_if(
      m_ElementsParams.begin(), m_ElementsParams.end(),
      [line_id](auto const &elt_ptr) { return elt_ptr->hasLine(line_id); });
  if (it != m_ElementsParams.end()) {
    elt_index = it - m_ElementsParams.begin();
    line_index = it->get()->getLineIndex(line_id);
  }
  return std::make_pair(elt_index, line_index);
}

std::pair<Int32, Int32>
CLMEltListVector::findElementIndex(const std::string &LineTagStr,
                                   CLine::EType linetype) const {
  Int32 elt_index = undefIdx;
  Int32 line_index = undefIdx;
  auto it = std::find_if(m_ElementsParams.begin(), m_ElementsParams.end(),
                         [&LineTagStr, &linetype](auto const &elt_ptr) {
                           return elt_ptr->hasLine(LineTagStr) &&
                                  (linetype == CLine::EType::nType_All ||
                                   elt_ptr->m_type == linetype);
                         });
  if (it != m_ElementsParams.end()) {
    elt_index = it - m_ElementsParams.begin();
    line_index = it->get()->getLineIndex(LineTagStr);
  }
  return std::make_pair(elt_index, line_index);
}

std::vector<std::pair<Int32, TInt32List>>
CLMEltListVector::getIgmLinesIndices() const {

  std::vector<std::pair<Int32, TInt32List>> indices;
  for (size_t elt_idx = 0; elt_idx < getNbElements(); ++elt_idx) {
    auto const &line_indices = m_ElementsParams[elt_idx]->m_asymLineIndices;
    if (!line_indices.empty())
      indices.push_back({elt_idx, line_indices});
  }
  return indices;
}

/**
 * \brief If outside lambda range, sets fitted amplitudes and errors to -1. If
 *inside, sets each line's fitted amplitude and error to -1 if line outside
 *lambda range, or amplitude to A * nominal amplitude and error to SNR * nominal
 *amplitude.
 **/
void CLMEltListVector::SetElementAmplitude(Int32 eltIndex, Float64 A,
                                           Float64 AStd) {
  m_ElementsParams[eltIndex]->setAmplitudes(A, AStd);
}

void CLMEltListVector::resetLambdaOffsets() {
  for (auto &ep : m_ElementsParams)
    ep->resetLambdaOffsets();
}

void CLMEltListVector::resetAmplitudeOffsets() {
  for (auto &ep : m_ElementsParams)
    ep->resetAmplitudeOffset();
}

void CLMEltListVector::resetElementsFittingParam(bool enableAmplitudeOffsets) {

  for (auto const &ep : m_ElementsParams) {
    ep->resetFittingParams();
    if (enableAmplitudeOffsets)
      ep->resetAmplitudeOffset();
  }
}

void CLMEltListVector::resetAsymfitParams() {
  for (auto const &ep : m_ElementsParams) {
    ep->resetAsymfitParams();
  }
}

void CLMEltListVector::setGlobalOutsideLambdaRangeFromSpectra() {
  for (size_t elt_idx = 0; elt_idx < getNbElements(); ++elt_idx) {
    m_ElementsParams[elt_idx]->m_globalOutsideLambdaRange =
        computeOutsideLambdaRange(elt_idx);
    for (size_t line_idx = 0; line_idx < m_ElementsParams[elt_idx]->size();
         ++line_idx)
      m_ElementsParams[elt_idx]->m_globalOutsideLambdaRangeList[line_idx] =
          computeOutsideLambdaRangeLine(elt_idx, line_idx);
  }
}

void CLMEltListVector::setAllAbsLinesFittable() {
  m_absLinesNoContinuum = false;
  for (auto const &elt_param_ptr : m_ElementsParams) {
    elt_param_ptr->m_absLinesNullContinuum = false; // reset all
  }
}

void CLMEltListVector::setAllAbsLinesNotFittable() {
  m_absLinesNoContinuum = true;
  for (auto const &elt_param_ptr : m_ElementsParams) {
    elt_param_ptr->m_absLinesNullContinuum = false; // reset all
    if (elt_param_ptr->GetElementType() == CLine::EType::nType_Absorption)
      elt_param_ptr->m_absLinesNullContinuum = true;
  }
}

void CLMEltListVector::setNullNominalAmplitudesNotFittable() {
  for (auto &elt_param_ptr : m_ElementsParams) {
    elt_param_ptr->m_nullNominalAmplitudes = false; // reset all
    if (!elt_param_ptr->m_globalOutsideLambdaRange)
      elt_param_ptr->setNullNominalAmplitudesNotFittable();
  }
}

void CLMEltListVector::setAbsLinesNullContinuumNotFittable(
    CSpcModelVectorPtr const &models) {
  for (size_t eIdx = 0; eIdx < getNbElements(); ++eIdx) {
    auto const &elt_param_ptr = m_ElementsParams[eIdx];
    elt_param_ptr->m_absLinesNullContinuum = false; // reset all
    if (elt_param_ptr->GetElementType() != CLine::EType::nType_Absorption ||
        elt_param_ptr->m_globalOutsideLambdaRange)
      continue;
    if (m_absLinesNoContinuum ||
        models->getMaxContinuumUnderElement(eIdx) <= 0.)
      elt_param_ptr->m_absLinesNullContinuum = true;
  }
}

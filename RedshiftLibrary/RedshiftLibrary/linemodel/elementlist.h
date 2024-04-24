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
#ifndef _REDSHIFT_LINEMODEL_ELEMENTLIST_
#define _REDSHIFT_LINEMODEL_ELEMENTLIST_

#include "RedshiftLibrary/linemodel/element.h"

namespace NSEpic {
class CLineModelSolution;
class CLineModelElementList {
private:
  std::vector<std::shared_ptr<CLineModelElement>> m_Elements;

public:
  TInt32List m_elementsDisabledIndexes;

  TInt32List GetModelValidElementsIndexes() const;
  TInt32List getValidElementIndices(CLine::EType lineTypeFilter) const;

  Float64 GetElementAmplitude(Int32 j) const;

  TInt32List getOverlappingElements(Int32 ind, const TInt32Set &excludedInd,
                                    Float64 redshift,
                                    Float64 overlapThres) const;

  Int32 GetModelValidElementsNDdl() const;

  std::vector<TInt32List> GetModelVelfitGroups(CLine::EType lineType) const;

  TInt32RangeList getlambdaIndexesUnderLines(
      const TInt32List &eIdx_list, const TInt32List &subeIdx_list,
      Float64 sigma_support, const CSpectrumSpectralAxis &spectralAxis,
      const TFloat64Range &lambdaRange, Float64 redshift) const;

  TInt32List getSupportIndexes(const TInt32List &EltsIdx) const;

  void addToSpectrumAmplitudeOffset(
      const CSpectrumSpectralAxis &spectralAxis,
      CSpectrumFluxAxis &modelfluxAxis, const TInt32List &eIdx_list = {},
      CLine::EType lineTypeFilter = CLine::EType::nType_All) const;

  bool IsElementIndexInDisabledList(Int32 index) const;
  void SetElementIndexesDisabledAuto();
  void ResetElementIndexesDisabled();

  Float64 getScaleMargCorrection(Int32 Eltidx = undefIdx) const;
  bool GetModelStrongEmissionLinePresent() const;
  bool GetModelHaStrongest() const;

  const TPolynomCoeffs &getPolynomCoeffs(Int32 eIdx) const;
  std::map<std::string, TInt32List>
  getFittingGroups(TInt32List EltsIdx = {},
                   CLine::EType lineTypeFilter = CLine::EType::nType_All) const;
  std::pair<Int32, Int32> findElementIndex(Int32 line_id) const;

  std::vector<std::shared_ptr<CLineModelElement>> &getElements() {
    return m_Elements;
  }

  const std::shared_ptr<const CLineModelElement> operator[](Int32 i) const {
    return m_Elements[i];
  }
  std::shared_ptr<CLineModelElement> &operator[](Int32 i) {
    return m_Elements[i];
  }

  // for range looping over elements
  std::vector<std::shared_ptr<CLineModelElement>>::iterator begin() {
    return m_Elements.begin();
  };
  std::vector<std::shared_ptr<CLineModelElement>>::iterator end() {
    return m_Elements.end();
  };
  std::vector<std::shared_ptr<CLineModelElement>>::const_iterator
  begin() const {
    return m_Elements.begin();
  };
  std::vector<std::shared_ptr<CLineModelElement>>::const_iterator end() const {
    return m_Elements.end();
  };

  Int32 size() const { return m_Elements.size(); }
  void push_back(const std::shared_ptr<CLineModelElement> &elt) {
    m_Elements.push_back(elt);
  }

  void debug(std::ostream &os) const {
    for (auto elt : m_Elements)
      elt->debug(os);
  }
};

class CLMEltListVector {
public:
  CLMEltListVector(CTLambdaRangePtrVector lambdaranges,
                   const std::shared_ptr<Int32> &curObs,
                   const CLineMap &restLineList,
                   ElementComposition element_composition);
  CLMEltListVector(CLineModelElementList eltlist,
                   const CLineMap &restLineList); // for unit test

  CLMEltListVector() = delete;

  std::pair<Int32, Int32> findElementIndex(Int32 line_id) const;
  std::pair<Int32, Int32>
  findElementIndex(const std::string &LineTagStr,
                   CLine::EType linetype = CLine::EType::nType_All) const;
  Int32 getNbElements() const { return m_ElementsParams.size(); }
  Float64 getScaleMargCorrection(Int32 Eltidx = undefIdx) const;

  TInt32List findElementTypeIndices(CLine::EType type) const;
  TInt32List getSupportIndexes(const TInt32List &EltsIdx) const;
  std::vector<std::pair<Int32, TInt32List>> getIgmLinesIndices() const;

  bool isOutsideLambdaRangeLine(Int32 elt_index, Int32 line_index) const {
    return m_globalOutsideLambdaRangeList[elt_index][line_index];
  }
  bool isOutsideLambdaRange(Int32 elt_index) const {
    return m_globalOutsideLambdaRange[elt_index];
  }
  const std::vector<bool> &getOutsideLambdaRangeList(Int32 elt_index) const;

  CLineModelElementList &getElementList() {
    return m_ElementsVector.at(*m_curObs);
  }
  const CLineModelElementList &getElementList() const {
    return m_ElementsVector.at(*m_curObs);
  }
  std::vector<TLineModelElementParam_ptr> &getElementParam() {
    return m_ElementsParams;
  }
  const std::vector<TLineModelElementParam_ptr> &getElementParam() const {
    return m_ElementsParams;
  }

  void SetElementAmplitude(Int32 eltIndex, Float64 A, Float64 AStd);

  void resetLambdaOffsets();
  void resetAmplitudeOffsets();
  void resetElementsFittingParam(bool enableAmplitudeOffsets);
  void resetAsymfitParams();

  void setGlobalOutsideLambdaRangeFromSpectra();

  Int32 GetModelNonZeroElementsNDdl() const;

private:
  std::vector<CLineModelElementList> m_ElementsVector;
  std::vector<TLineModelElementParam_ptr> m_ElementsParams;
  CTLambdaRangePtrVector m_lambdaRanges;
  Int32 m_nbObs;
  std::shared_ptr<Int32> m_curObs;
  const CLineMap &m_RestLineList;

  std::vector<bool> m_globalOutsideLambdaRange;
  std::vector<std::vector<bool>> m_globalOutsideLambdaRangeList;

  void AddElementParam(CLineVector lines);
  void fillElements();
  void LoadCatalog();
  void LoadCatalogOneLineByElement();
  void LoadCatalogOneMultiline();
  void LoadCatalogTwoMultilinesAE();
  bool computeOutsideLambdaRangeLine(Int32 elt_index, Int32 line_index);
  bool computeOutsideLambdaRange(Int32 elt_index);
};

} // namespace NSEpic
#endif

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
#include <cstdarg>
#include <iostream>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/ruleStrongHigherThanWeak.h"
#include "RedshiftLibrary/log/log.h"
using namespace NSEpic;
using namespace std;

CRuleStrongHigherThanWeak::CRuleStrongHigherThanWeak()
    : m_LineType(CLine::EType::nType_All) {}

void CRuleStrongHigherThanWeak::SetUp(bool EnabledArgument, ...) {
  Name = "strongWeak";
  Enabled = EnabledArgument;
  va_list Arguments;
  va_start(Arguments, EnabledArgument);
  m_LineType = va_arg(Arguments, CLine::EType);
  va_end(Arguments);
}

/** \brief Verify that "stronger lines have higher amplitudes than weaker lines"
 *rule is applicable, and then apply it. If the maximum amplitude for strong
 *lines of the specified type is -1, do nothing. For each weak line of the
 *specified type: Find the first element that contains that weak line Get the
 *index of the entry in that element that corresponds to the weak line If the
 *indexed entry IsOutsideLambdaRange, go for the next weak line Get the
 *parameters for the entry Limit the amplitude of the entry to the maximum
 *amplitude for strong lines
 **/
void CRuleStrongHigherThanWeak::Correct(
    CLMEltListVector &LineModelElementList) {
  auto const [minStrongEltIndex, minStrongLineIndex] =
      FindLowestStrongLineIndex(LineModelElementList);

  // If no strong line was found, there is no correction to apply
  if (minStrongEltIndex == undefIdx) {
    return;
  }

  // Access strong line infos
  auto &minElement_ptr =
      LineModelElementList.getElementParam()[minStrongEltIndex];
  Float64 erStrong = minElement_ptr->GetFittedAmplitudeStd(minStrongLineIndex);
  Float64 ampStrong = minElement_ptr->GetFittedAmplitude(minStrongLineIndex);
  std::string nameStrong =
      minElement_ptr->GetLines()[minStrongLineIndex].GetName();
  Float64 maxAmp = maxAmplitude(ampStrong, erStrong);

  for (Int32 iElement = 0;
       iElement < LineModelElementList.getElementList().size(); iElement++) {
    auto &element_ptr = LineModelElementList.getElementList()[iElement];
    auto const &element_param_ptr = element_ptr->getElementParam();

    // Consider only desired line types
    if (element_param_ptr->GetElementType() != m_LineType)
      continue;
    correctLineModelElement(*element_ptr, maxAmp, nameStrong);
  }
  Log.LogInfo(Logs);
}

void CRuleStrongHigherThanWeak::correctLineModelElement(
    CLineModelElement &element, Float64 maxAmplitude,
    const std::string &nameStrong) {

  const auto &elt_param_ptr = element.getElementParam();
  for (Int32 iLine = 0; iLine != element.GetSize(); ++iLine) {
    if (element.IsOutsideLambdaRange(iLine))
      continue;
    auto const &line = elt_param_ptr->GetLines()[iLine];
    if (line.IsWeak()) {
      Float64 fittedAmplitude = elt_param_ptr->GetFittedAmplitude(iLine);
      bool limited = elt_param_ptr->LimitFittedAmplitude(iLine, maxAmplitude);
      if (limited)
        constructLogMsg(line.GetName(), nameStrong, fittedAmplitude,
                        maxAmplitude);
    }
  }
}

Float64 CRuleStrongHigherThanWeak::maxAmplitude(Float64 ampStrong,
                                                Float64 erStrong) {
  // Calculates maximal amplitude, taking noise into account. Some possible
  // methods are : Method 0 : no noise taken into acccount Float64 maxAmplitude
  // = (coeff*ampStrong);

  // Method 1 : using Strong line noise to be more tolerant

  // Method 2 : using Strong line noise and Weak line noise to correct with
  // a ratio
  //      Float64 maxAmplitude = ampStrong; //default value
  //      if(erB>0.0 && erB<erA && erA>0.0)
  //      {
  //          maxAmplitude = (coeff*ampStrong)*(erA/erB);
  //      }else{
  //          maxAmplitude = (coeff*ampStrong);
  //      }
  //

  Float64 COEFF = 1.0;
  Float64 N_SIGMA = 1.0;
  if (isnan(erStrong)) {
    THROWG(INTERNAL_ERROR, Formatter() << "CRuleStrongHigherThanWeak"
                                       << __func__ << " erStrong is Nan");
  }
  if (isnan(ampStrong)) {
    THROWG(INTERNAL_ERROR, Formatter() << "CRuleStrongHigherThanWeak"
                                       << __func__ << " ampStrong is Nan");
  }

  // Method 1 is used here
  return (COEFF * ampStrong) + COEFF * (erStrong * N_SIGMA);
}

void CRuleStrongHigherThanWeak::constructLogMsg(const std::string &nameWeak,
                                                const std::string &strongName,
                                                Float64 fittedAmplitude,
                                                Float64 maxAmplitude) {
  if (Logs.size() == 0) {
    std::string strTmp0 =
        boost::str((boost::format("correct - %-10s") % "STRONG_WEAK"));
    Logs.append(strTmp0.c_str());
  }
  std::string strTmp =
      boost::str((boost::format("\n\tlineWeak=%-10s, lineStrong=%-10s, "
                                "previousAmp=%.4e, correctedAmp=%.4e") %
                  nameWeak % strongName % fittedAmplitude % maxAmplitude));
  Logs.append(strTmp.c_str());
}

bool CRuleStrongHigherThanWeak::Check(CLMEltListVector &LineModelElementList) {
  return false;
}

/**
 * \brief Returns the maximum amplitude between superstrong lines within the
 *support of m_Elements. The referenced er argument will hold the error sigma
 *for the same element.
 **/
Float64 CRuleStrongHigherThanWeak::FindHighestStrongLineAmp(
    Int32 linetype, Float64 &er, std::string &name,
    CLMEltListVector &LineModelElementList) {
  Float64 maxi = -1.0;
  for (Int32 iedx = 0; iedx < LineModelElementList.getElementList().size();
       iedx++) {
    const auto &element_ptr = LineModelElementList.getElementList()[iedx];
    const auto &element_param_ptr =
        LineModelElementList.getElementParam()[iedx];
    for (Int32 iLineStrong = 0; iLineStrong != element_ptr->GetSize();
         ++iLineStrong) {
      auto const &lineStrong = element_param_ptr->GetLines()[iLineStrong];
      if (element_ptr->IsOutsideLambdaRange(iLineStrong) ||
          lineStrong.GetForce() != CLine::EForce::nForce_Strong ||
          element_param_ptr->GetElementType() != m_LineType) {
        continue;
      }

      Float64 ampStrong =
          element_ptr->getElementParam()->GetFittedAmplitude(iLineStrong);
      Float64 erStrong =
          element_ptr->getElementParam()->GetFittedAmplitudeStd(iLineStrong);
      if (maxi < ampStrong /*&& lineSnr>validSNRCut*/) {
        maxi = ampStrong;
        er = erStrong;
        name = lineStrong.GetName();
      }
    }
  }
  return maxi;
}

std::pair<Int32, Int32> CRuleStrongHigherThanWeak::FindLowestStrongLineIndex(
    const CLMEltListVector &LineModelElementList) {
  Float64 amplitudeMin = INFINITY;
  Int32 iElementMin = undefIdx;
  Int32 iLineMin = undefIdx;
  for (Int32 iElement = 0;
       iElement < LineModelElementList.getElementList().size(); iElement++) {
    const auto &element_ptr = LineModelElementList.getElementList()[iElement];
    const auto &element_param_ptr =
        LineModelElementList.getElementParam()[iElement];
    if (element_param_ptr->GetElementType() != m_LineType)
      continue;

    for (Int32 iLine = 0; iLine != element_ptr->GetSize(); ++iLine) {
      auto const &line = element_param_ptr->GetLines()[iLine];
      if (element_ptr->IsOutsideLambdaRange(iLine))
        continue;

      if (line.IsStrong()) {
        Float64 lineAmplitude =
            element_ptr->getElementParam()->GetFittedAmplitude(iLine);
        if (lineAmplitude < amplitudeMin) {
          amplitudeMin = lineAmplitude;
          iLineMin = iLine;
          iElementMin = iElement;
        }
      }
    }
  }
  return std::make_pair(iElementMin, iLineMin);
}

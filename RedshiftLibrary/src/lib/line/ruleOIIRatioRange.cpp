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
#include "RedshiftLibrary/line/ruleOIIRatioRange.h"
#include "RedshiftLibrary/log/log.h"
#include <cstdarg>
#include <iostream>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace NSEpic;
using namespace std;

CRuleRatioRange::CRuleRatioRange()
    : m_LineType(0), m_LineA(""), m_LineB(""), m_Coefficient(0) {}

void CRuleRatioRange::SetUp(bool EnabledArgument, ...) {
  Name = "ratiorange";
  Enabled = EnabledArgument;
  va_list Arguments;
  va_start(Arguments, EnabledArgument);
  m_LineType = va_arg(Arguments, Int32);
  m_LineA = std::string(va_arg(Arguments, const char *));
  m_LineB = std::string(va_arg(Arguments, const char *));
  m_Coefficient = va_arg(Arguments, Float64);
  va_end(Arguments);
}

bool CRuleRatioRange::Check(CLineModelElementList &LineModelElementList) {
  return false;
}

/**
 * For two distinct lines, if neither IsOutsideLambdaRange, and their amplitudes
 *are beyond a range (considering coeff), SetFittedAmplitude of each with
 *corrected values.
 **/
void CRuleRatioRange::Correct(CLineModelElementList &LineModelElementList) {
  Int32 iA = LineModelElementList.findElementIndex(m_LineA, m_LineType);
  if (iA == -1 || LineModelElementList[iA]->GetSize() < 1) {
    Log.LogDebug("Rule %s: line %s not found or line has size <1", Name.c_str(),
                 m_LineA.c_str());
    return;
  }

  Int32 iB = LineModelElementList.findElementIndex(m_LineB, m_LineType);
  if (iB == -1 || LineModelElementList[iB]->GetSize() < 1) {
    Log.LogDebug("Rule %s: line %s not found or line has size <1", Name.c_str(),
                 m_LineB.c_str());
    return;
  }

  if (iA == iB)
    return;

  if (LineModelElementList[iA]->IsOutsideLambdaRange() ||
      LineModelElementList[iB]->IsOutsideLambdaRange())
    return;

  Float64 ampA = LineModelElementList[iA]->GetFittedAmplitude(0);
  Float64 erA = LineModelElementList[iA]->GetFittedAmplitudeErrorSigma(0);
  Float64 ampB = LineModelElementList[iB]->GetFittedAmplitude(0);
  Float64 erB = LineModelElementList[iB]->GetFittedAmplitudeErrorSigma(0);
  Int32 i1 = iA;
  Int32 i2 = iB;
  Float64 amp1 = ampA;
  Float64 er1 = erA;
  Float64 amp2 = ampB;
  Float64 er2 = erB;
  if (std::abs(ampA) > std::abs(ampB * m_Coefficient)) {
    i1 = iA;
    i2 = iB;
    amp1 = ampA;
    er1 = erA;
    amp2 = ampB;
    er2 = erB;
  } else {
    if (std::abs(ampB) > std::abs(ampA * m_Coefficient)) {
      i1 = iB;
      i2 = iA;
      amp1 = ampB;
      er1 = erB;
      amp2 = ampA;
      er2 = erA;
    } else {
      return;
    }
  }
  Float64 R = m_Coefficient;
  Float64 w1 = 0.0;
  if (er1 != 0.0)
    w1 = 1.0 / (er1 * er1);

  Float64 w2 = 0.0;
  if (er2 != 0.0)
    w2 = 1.0 / (er2 * er2 * R * R);

  Float64 corrected1 = (amp1 * w1 + amp2 * w2 * R) / (w1 + w2);
  Float64 corrected2 = corrected1 / R;

  // log the correction
  Float64 correctedA = corrected1;
  Float64 correctedB = corrected2;
  if (i1 == iB) {
    correctedA = corrected2;
    correctedB = corrected1;
  }
  constructLogMsg(m_LineA, ampA, correctedA);
  constructLogMsg(m_LineB, ampB, correctedB);

  LineModelElementList[i1]->SetElementAmplitude(corrected1, er1);
  LineModelElementList[i2]->SetElementAmplitude(corrected2, er2);
}

void CRuleRatioRange::constructLogMsg(const std::string &lineStrA, Float64 ampA,
                                      Float64 correctedA) {

  if (Logs.size() == 0) {
    std::string strTmp0 =
        boost::str((boost::format("correct - %-10s") % "RATIO_RANGE"));
    Logs.append(strTmp0.c_str());
  }
  std::string strTmp1 = boost::str(
      (boost::format("\n\tline=%-10s, previousAmp=%.4e, correctedAmp=%.4e") %
       lineStrA % ampA % correctedA));
  Logs.append(strTmp1.c_str());
}

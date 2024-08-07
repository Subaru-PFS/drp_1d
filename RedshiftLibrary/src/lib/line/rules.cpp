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
#include <fstream>

#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/line/rules.h"
#include "RedshiftLibrary/operator/linedetection.h"

using namespace NSEpic;
using namespace std;

CRules::CRules(CSpectrum &spc, CLineDetectedCatalog &detectedCatalog,
               CLineCatalog &restCatalog, TFloat64Range &lambdaRange,
               Float64 winsize)
    : m_spc(spc), m_DetectedCatalog(detectedCatalog),
      m_RestCatalog(restCatalog), m_lambdaRange(lambdaRange),
      m_winsize(winsize) {}

Int32 CRules::check(
    Float64 z, CLineMatchingResult::TSolutionSet &matchingSolutionSet) const {
  Int32 intval = 0;

  bool bval;
  bval = checkRule01(z, matchingSolutionSet);
  if (!bval)
    intval += 1;

  bval = checkRule02(z, matchingSolutionSet);
  if (!bval)
    intval += 10;

  bval = checkRule03(z, matchingSolutionSet);
  if (!bval)
    intval += 100;

  return intval;
}

/**
 * This rule applies when only weak rest lines have been identified in the
 * matching solution. In this case, the absence of strong rest lines must be
 * explained by:
 *  - the strong lines being out of the wavelength range_error
 *  - the string lines being unindentified due to high noise
 *
 */
bool CRules::checkRule01(
    Float64 z, CLineMatchingResult::TSolutionSet &matchingSolutionSet) const {

  // check if only weak lines are in this solution set
  for (const auto &mSolution : matchingSolutionSet)
    if (mSolution.RestLine.IsStrong())
      return true;

  // check if the absence of strong lines is justified by the wavelength range
  CLineMap const strongRestLineList = m_RestCatalog.GetFilteredList(
      CLine::EType::nType_Emission, CLine::EForce::nForce_Strong);
  CLineMap strongRestLinesInsideLambdaRangeList;

  for (const auto &[id, line] : strongRestLineList) {
    Float64 lambda = line.GetPosition() * (1 + z);
    if (isLineInsideRange(lambda))
      strongRestLinesInsideLambdaRangeList[id] = line;
  }

  if (strongRestLinesInsideLambdaRangeList.empty())
    return true;

  // check if the absence of the remaining strong lines is explained by noise
  CLineDetection lineDetection;
  // estimate the weak lines SNR
  Float64 maxSnrWeak = 0.0;
  Float64 maxNoiseWeak = 0.0;

  for (const auto &mSolution : matchingSolutionSet) {
    Float64 const lambda = mSolution.RestLine.GetPosition() * (1 + z);
    TFloat64Range const lambdarange(lambda - m_winsize / 2.0,
                                    lambda + m_winsize / 2.0);
    TInt32Range const range =
        m_spc.GetSpectralAxis().GetIndexesAtWaveLengthRange(lambdarange);

    Float64 flux = 0.0;
    Float64 noiseWeak = 0.0;
    Float64 snrWeak = lineDetection.ComputeFluxes(
        m_spc, m_winsize, range, TFloat64List(), &flux, &noiseWeak);

    maxSnrWeak = std::max(maxSnrWeak, snrWeak);
    maxNoiseWeak = std::max(maxNoiseWeak, noiseWeak);
  }
  // estimate the strong lines SNR
  for (const auto &[_, line] : strongRestLineList) {
    Float64 const lambda = line.GetPosition() * (1 + z);
    TFloat64Range const lambdarange(lambda - m_winsize / 2.0,
                                    lambda + m_winsize / 2.0);
    TInt32Range const range =
        m_spc.GetSpectralAxis().GetIndexesAtWaveLengthRange(lambdarange);

    Float64 flux = 0.0;
    Float64 noise = 0.0;
    Float64 snrStrong = lineDetection.ComputeFluxes(
        m_spc, m_winsize, range, TFloat64List(), &flux, &noise);

    if (snrStrong < maxSnrWeak && noise < maxNoiseWeak)
      return false;
  }
  return true;
}

/**
 * This rule applies when one of OIII lines have been detected.
 *
 */
bool CRules::checkRule02(
    Float64 z, CLineMatchingResult::TSolutionSet &matchingSolutionSet) const {

  // check if the OIII doublet is in this solution set
  Int32 founda = 0;
  Int32 foundb = 0;

  for (const auto &mSolution : matchingSolutionSet) {
    std::string const &name = mSolution.RestLine.GetName();
    std::size_t const foundstra = name.find(linetags::oIIIa_em);
    if (foundstra != std::string::npos)
      founda++;

    std::size_t const foundstrb = name.find(linetags::oIIIb_em);
    if (foundstrb != std::string::npos) {
      // check if OIIIa would be in the wavelength range
      Float64 const lambda = getRestLineLambda(linetags::oIIIa_em) * (1 + z);
      if (isLineInsideRange(lambda, true))
        foundb++;
    }
  }
  if (foundb == 1 && founda == 0) {
    return false;
  }

  return true;
}

/**
 * This rule applies when Hbeta has been detected.
 *
 */
bool CRules::checkRule03(
    Float64 z, CLineMatchingResult::TSolutionSet &matchingSolutionSet) const {

  // check if the Hbeta doublet is in this solution set
  Int32 foundHbeta = 0;
  Int32 foundHalpha = 0;

  // Float64 z =
  // CLineMatchingResult::GetMeanRedshiftSolution(matchingSolutionSet);
  for (const auto &mSolution : matchingSolutionSet) {
    std::string const &name = mSolution.RestLine.GetName();
    std::size_t foundstr = name.find(linetags::hbeta_em);
    if (foundstr != std::string::npos) {
      // check if Halpha would be in the wavelength range
      Float64 const lambda = getRestLineLambda(linetags::halpha_em) * (1 + z);
      if (isLineInsideRange(lambda, true))
        foundHbeta++;
    }
    foundstr = name.find(linetags::halpha_em);
    if (foundstr != std::string::npos)
      foundHalpha++;
  }
  if (foundHbeta == 1 && foundHalpha != 1)
    return false;

  return true;
}

Float64 CRules::getRestLineLambda(std::string nametag) const {
  CLineMap const restLineList = m_RestCatalog.GetFilteredList();
  for (auto const &[id, line] : restLineList) {
    std::string const &name = line.GetName();
    std::size_t const foundstr = name.find(nametag.c_str());
    if (foundstr != std::string::npos)
      return line.GetPosition();
  }
  return NAN;
}

bool CRules::isLineInsideRange(Float64 lambda, bool winsize) const {
  Float64 start = m_lambdaRange.GetBegin() + winsize ? m_winsize : 0.;
  Float64 end = m_lambdaRange.GetEnd() + winsize ? m_winsize : 0.;
  return (lambda >= start && lambda <= end);
}

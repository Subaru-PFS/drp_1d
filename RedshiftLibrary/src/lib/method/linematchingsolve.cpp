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
#include "RedshiftLibrary/method/linematchingsolve.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/extremum/extremum.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/method/linematchingsolveresult.h"
#include "RedshiftLibrary/operator/linedetection.h"
#include "RedshiftLibrary/operator/linedetectionresult.h"
#include "RedshiftLibrary/operator/linematching.h"
#include "RedshiftLibrary/operator/linematchingresult.h"
#include "RedshiftLibrary/operator/peakdetection.h"
#include "RedshiftLibrary/operator/peakdetectionresult.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"

using namespace NSEpic;
using namespace std;

/**
 * \brief This constructor will attribute values to this method's parameters
 * with default values.
 */
CLineMatchingSolve::CLineMatchingSolve() : CObjectSolve("LineMatchingSolve") {
  Log.LogDebug("CLineMatchingSolve::CLineMatchingSolve()");

  // Peak Detection
  m_winsize = 250.0;
  m_cut = 5.0;
  m_detectioncut = 5.0;
  m_detectionnoiseoffset = 0.0; /**\var Useful for simulation data. */
  m_strongcut = 2.0;
  m_minsize = 3;
  m_maxsize = 70;
  m_enlargeRate = 2.0;

  // Line Matching
  m_disablegaussianfitqualitycheck = false;
  m_minMatchNum = 1;
  m_tol = 0.002;
}

/**
 * Empty destructor.
 */
CLineMatchingSolve::~CLineMatchingSolve() {
  Log.LogDebug("CLineMatchingSolve::~CLineMatchingSolve()");
}

/**
 * \brief Attempts to find a matching line using the input parameters.
 * This algorithm will create and run a CPeakDetection object, using the
 * parameters as input.
 *
 * The input spectrum is copied to the local scope.
 *
 * If the lines' type are absorption lines, the spectrum has its fluxes
 * inverted.
 *
 * This local spectrum is passed to a lineDetection object.
 *
 * If the dynamic option is false, the following algorithm will be run:
 *
 * - If peaks are found, it will try to find lines in the peaks found, using a
 * CLineDetection object (please note that this class is defined in
 * linedetectionresult.h).
 *
 * - If lines are found, it will try to match these lines to the catalogue of
 * known lines and determine the redshift using this information.
 *
 * If the dynamic option is true, the line detection is run several times, with
 * the parameters being varied withing physically meaningful values. When either
 * a threshold number of peaks is detected, or all parameters are exhaustively
 * searched, the algorithm continues as normal.
 */
std::shared_ptr<CSolveResult> CLineMatchingSolve::compute() {

  auto &resultStore = *Context.GetResultStore();
  auto const &paramStore = Context.GetParameterStore();
  auto const &spc = *Context.GetInputContext()->GetSpectrum();
  auto const &restLineCatalog =
      *Context.GetInputContext()->GetLineCatalog(m_category, m_name);

  auto lineType = CLine::EType::nType_Emission;
  std::string linetypeStr = "E";

  Log.LogDebug("Attempting to load parameters from parameter JSON.");
  {
    m_cut = paramStore->GetScoped<Float64>("linematching.cut");
    m_detectioncut =
        paramStore->GetScoped<Float64>("linematching.detectioncut");
    m_detectionnoiseoffset =
        paramStore->GetScoped<Float64>("linematching.detectionnoiseoffset");
    m_disablegaussianfitqualitycheck = paramStore->GetScoped<bool>(
        "linematching.disablegaussianfitqualitycheck");
    m_dynamicLinematching =
        paramStore->GetScoped<bool>("linematching.dynamicLinematching");
    m_enlargeRate = paramStore->GetScoped<Float64>("linematching.enlargeRate");
    linetypeStr = paramStore->GetScoped<std::string>("linematching.linetype");
    m_maxsize = paramStore->GetScoped<Float64>("linematching.maxsize");
    m_minMatchNum = paramStore->GetScoped<Int32>("linematching.minMatchNum");
    m_minsize = paramStore->GetScoped<Float64>("linematching.minsize");
    m_strongcut = paramStore->GetScoped<Float64>("linematching.strongcut");
    m_tol = paramStore->GetScoped<Float64>("linematching.tol");
    m_winsize = paramStore->GetScoped<Float64>("linematching.winsize");
    lineType = CLine::EType::nType_All;
    if (linetypeStr == "Absorption") {
      lineType = CLine::EType::nType_Absorption;
    }
    if (linetypeStr == "Emission") {
      lineType = CLine::EType::nType_Emission;
    }
  }

  Log.LogDebug("Final parameters read:");
  {
    Log.LogDebug(Formatter() << "m_cut =" << m_cut);
    if (m_dynamicLinematching) {
      Log.LogDebug(
          "m_dynamicLinematching is true (any value, except 0, on json)");
    } else {
      Log.LogDebug("m_dynamicLinematching is false (value 0 on json)");
    }
    Log.LogDebug(Formatter() << "m_detectioncut =" << m_detectioncut);
    Log.LogDebug(Formatter()
                 << "m_detectionnoiseoffset =" << m_detectionnoiseoffset);
    if (m_disablegaussianfitqualitycheck) {
      Log.LogDebug("m_disablegaussianfitqualitycheck is true (any value, "
                   "except 0, in json)");
    } else {
      Log.LogDebug(
          "m_disablegaussianfitqualitycheck is false (value 0 in json)");
    }
    Log.LogDebug(Formatter() << "m_enlargeRate = " << m_enlargeRate);
    Log.LogDebug(Formatter() << "linetype = " << linetypeStr);
    Log.LogDebug(Formatter() << "m_minMatchNum = " << m_minMatchNum);
    Log.LogDebug(Formatter() << "m_minsize = " << m_minsize);
    Log.LogDebug(Formatter() << "m_maxsize = " << m_maxsize);
    Log.LogDebug(Formatter() << "m_strongcut = " << m_strongcut);
    Log.LogDebug(Formatter() << "m_tol = " << m_tol);
    Log.LogDebug(Formatter() << "m_winsize = " << m_winsize);
  }

  CPeakDetection peakDetection(m_winsize, m_detectioncut, 1, m_enlargeRate,
                               m_detectionnoiseoffset);
  CSpectrum _spc = spc;
  if (lineType == CLine::EType::nType_Absorption) {
    Log.LogDebug("lineType == CLine::EType::nType_Absorption");
    _spc.InvertFlux();
  }

  if (lineType == CLine::EType::nType_Emission) {
    Log.LogDebug("lineType == CLine::EType::nType_Emission");
  }

  auto peakDetectionResult = peakDetection.Compute(_spc);
  if (peakDetectionResult) {
    Log.LogDebug(Formatter()
                 << "Storing " << peakDetectionResult->PeakList.size()
                 << " peaks from PeakList in the result store.");
    resultStore.StoreScopedGlobalResult("peakdetection", peakDetectionResult);
  } else {
    THROWG(ErrorCode::INTERNAL_ERROR, "No peak detected");
  }

  // Since we detected at least one peak, try to detect lines related to those
  // peaks.
  Int32 bestMatchingNumber = -1;
  Float64 bestRedshift = -1.0;
  Int32 currentNumberOfPeaks = -1;
  Float64 cutBest = m_cut;
  auto cutCurrent = m_cut;
  Float64 cutMaximum = 1e1;
  Float64 cutMinimum = 1e-1;
  Float64 cutStep = (cutMaximum - cutMinimum) / 5;
  Float64 minimumFwhhBest = m_minsize;
  Float64 minimumFwhhCurrent = m_minsize;
  Float64 minimumFwhhMinimum = 1e-3;
  Float64 minimumFwhhMaximum = 1e1;
  Float64 minimumFwhhStep = (minimumFwhhMaximum - minimumFwhhMinimum) / 5;
  Int32 minimumNumberOfPeaks = 200;
  bool newValues = false;
  Int32 numberOfPeaksBest = -1;
  bool numberOfPeaksBestBypass = false;
  Int32 previousNumberOfPeaks = -1;
  Float64 strongcutBest = m_strongcut;
  auto strongcutCurrent = m_strongcut;
  Float64 strongcutMaximum = 1e1;
  Float64 strongcutMinimum = 1e-2;
  Float64 strongcutStep = (strongcutMaximum - strongcutMinimum) / 5;
  Float64 winsizeBest = m_winsize;
  auto winsizeCurrent = m_winsize;
  Float64 winsizeMaximum = 1e4;
  Float64 winsizeMinimum = 1e3;
  Float64 winsizeStep = (winsizeMaximum - winsizeMinimum) / 5;

  Int32 cmptMax = 100;
  Int32 iCmpt = 0;
  while (iCmpt < cmptMax) {
    Log.LogDebug("*************************");
    Log.LogDebug(Formatter() << "cutCurrent == " << cutCurrent);
    Log.LogDebug(Formatter() << "minimumFwhhCurrent == " << minimumFwhhCurrent);
    Log.LogDebug(Formatter() << "strongcutCurrent == " << strongcutCurrent);
    Log.LogDebug(Formatter() << "winsizeCurrent == " << winsizeCurrent);
    CLineDetection lineDetection(lineType, cutCurrent, strongcutCurrent,
                                 winsizeCurrent, minimumFwhhCurrent, m_maxsize,
                                 m_disablegaussianfitqualitycheck);
    auto lineDetectionResult = lineDetection.Compute(
        _spc, m_lambdaRange, peakDetectionResult->PeakList,
        peakDetectionResult->EnlargedPeakList);
    if (lineDetectionResult || !m_dynamicLinematching) {
      currentNumberOfPeaks = lineDetectionResult->LineCatalog.GetList().size();
      Log.LogDebug(Formatter()
                   << "Found " << currentNumberOfPeaks << " peaks.");
      if (currentNumberOfPeaks >= minimumNumberOfPeaks ||
          numberOfPeaksBestBypass || !m_dynamicLinematching) {
        Log.LogDebug(Formatter()
                     << "Storing "
                     << lineDetectionResult->LineCatalog.GetList().size()
                     << " lines from lineDetection in the result store.");
        resultStore.StoreScopedGlobalResult("lineCatalog", lineDetectionResult);
        // Since we know at least one peak that corresponds to a line, let's try
        // to match to a catalogued template.
        CLineMatching lineMatching;
        Log.LogDebug("Now starting linematching");
        auto lineMatchingResult = lineMatching.Compute(
            lineDetectionResult->LineCatalog, restLineCatalog, m_redshifts,
            m_minMatchNum, m_tol, lineType);
        if (lineMatchingResult) {
          lineMatchingResult->FilterWithRules(_spc, m_lambdaRange, m_winsize);
          Log.LogDebug(Formatter()
                       << "CLineMatching yielded "
                       << lineMatchingResult->SolutionSetList.size()
                       << " sets of solutions and "
                       << lineMatchingResult->FilteredSolutionSetList.size()
                       << " sets of filtered solutions.");
          // Store matching results
          resultStore.StoreScopedGlobalResult("linematching",
                                              lineMatchingResult);
          lineMatchingResult->GetBestRedshift(bestRedshift, bestMatchingNumber);
          Log.LogDebug(Formatter()
                       << "bestRedshift == " << bestRedshift
                       << ", bestMatchingNumber == " << bestMatchingNumber);
          if (bestRedshift != -1.0) {
            Log.LogDebug(
                "return std::shared_ptr<const CLineMatchingSolveResult>( new "
                "CLineMatchingSolveResult() );");
            return std::make_shared<CLineMatchingSolveResult>();
          } else if (!m_dynamicLinematching) {
            return std::make_shared<CLineMatchingSolveResult>();
          }
        } // lineMatchingResult
        else if (!m_dynamicLinematching) {
          return std::make_shared<CLineMatchingSolveResult>();
        }
      } // minimumNumberOfPeaks
    }   // lineDetectionResult
    // update parameters
    {
      if (currentNumberOfPeaks >= numberOfPeaksBest) {
        numberOfPeaksBest = currentNumberOfPeaks;
        cutBest = cutCurrent;
        minimumFwhhBest = minimumFwhhCurrent;
        strongcutBest = strongcutCurrent;
        winsizeBest = winsizeCurrent;
      }
      if (currentNumberOfPeaks < previousNumberOfPeaks && !newValues)
        THROWG(ErrorCode::INTERNAL_ERROR,
               "Dynamic logic failed - number of peaks is falling.");

      previousNumberOfPeaks = currentNumberOfPeaks;
      if (!m_dynamicLinematching)
        THROWG(ErrorCode::INTERNAL_ERROR, "No result found");

      newValues = false;
      if (minimumFwhhCurrent > minimumFwhhMinimum) {
        minimumFwhhCurrent =
            std::max(minimumFwhhMinimum, minimumFwhhCurrent - minimumFwhhStep);
        continue;
      }
      newValues = true;
      minimumFwhhCurrent = minimumFwhhMaximum;
      if (cutCurrent >= cutMinimum) {
        cutCurrent -= cutStep;
        continue;
      }
      cutCurrent = strongcutCurrent;
      if (strongcutCurrent > strongcutMinimum) {
        strongcutCurrent =
            std::max(strongcutMinimum, strongcutCurrent - strongcutStep);
        continue;
      }
      strongcutCurrent = strongcutMaximum;
      if (winsizeCurrent <= winsizeMaximum) {
        winsizeCurrent += winsizeStep;
        continue;
      }
      winsizeCurrent = winsizeMinimum;
      // if all checks failed, dynamic system failed.
      Log.LogDebug(
          "Dynamic logic: all parameters searched. Using best values.");
      {
        cutCurrent = cutBest;
        minimumFwhhCurrent = minimumFwhhBest;
        strongcutCurrent = strongcutBest;
        winsizeCurrent = winsizeBest;
        numberOfPeaksBestBypass = true;
        continue;
      }
    }
    Log.LogDebug("Dynamically retrying linematching.");
    iCmpt++;
  }; // while

  if (iCmpt == cmptMax) {
    Flag.warning(WarningCode::LINEMATCHING_REACHED_ENDLOOP,
                 Formatter()
                     << "CMethodLineMatchingSolve::" << __func__
                     << ": Stopped the linematching dynamic cut loop...");
  }
  return std::make_shared<CLineMatchingSolveResult>();
}

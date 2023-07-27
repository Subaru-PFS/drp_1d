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
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/log/log.h"

#include <algorithm> // std::sort
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;

/**
 * @brief CLineCatalogsTplRatio::GetRestLinesList
 * @param index
 * WARNING: ismCoeff not applied on the restlines provided by that function.
 */
CLineMap CLineCatalogsTplRatio::GetRestLinesList(Int32 index) const {
  auto const typeFilter = CLine::EType::nType_All;
  auto const forceFilter = CLine::EForce::nForce_All;

  CLineMap restLineList =
      m_lineRatioCatalogs[index].GetFilteredList(typeFilter, forceFilter);
  return restLineList;
}

TFloat64List CLineCatalogsTplRatio::getCatalogsPriors() const {
  TFloat64List catalogsPriors;
  catalogsPriors.reserve(m_lineRatioCatalogs.size());
  for (CLineRatioCatalog cat : m_lineRatioCatalogs)
    catalogsPriors.push_back(cat.getPrior());
  return catalogsPriors;
}

bool CLineCatalogsTplRatio::GetCatalogVelocities(Int32 idx, Float64 &elv,
                                                 Float64 &alv) const {
  // TODO generic velocity groups : there should not be hardcoded values, this
  // should return a map
  elv = m_lineRatioCatalogs[idx].getVelocity("em_vel");
  alv = m_lineRatioCatalogs[idx].getVelocity("abs_vel");
  return true;
}

std::vector<TFloat64Map> CLineCatalogsTplRatio::InitLineCorrespondingAmplitudes(
    Int32 enableISMCalzetti,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti) const {

  std::vector<TFloat64Map> lineCatalogCorrespondingNominalAmp(
      GetCatalogsCount(), TFloat64Map());

  for (Int32 iCatalog = 0; iCatalog != GetCatalogsCount(); ++iCatalog) {
    auto const &catalog = GetCatalog(iCatalog);
    for (auto const &[id, line] : catalog.GetList()) {
      Float64 nominalAmp = line.GetNominalAmplitude();
      Float64 const restLambda = line.GetPosition();
      if (enableISMCalzetti) {
        Float64 dustCoeff = ismCorrectionCalzetti->GetDustCoeff(
            catalog.getIsmIndex(), restLambda);
        nominalAmp *= dustCoeff;
      }
      lineCatalogCorrespondingNominalAmp[iCatalog][id] = nominalAmp;
    }
  }

  // Now log the linesCorrespondingNominalAmp
  logLineNominalAmp(enableISMCalzetti, lineCatalogCorrespondingNominalAmp,
                    ismCorrectionCalzetti);

  return lineCatalogCorrespondingNominalAmp;
}

// only logging
void CLineCatalogsTplRatio::logLineNominalAmp(
    bool enableISMCalzetti,
    const std::vector<std::map<Int32, Float64>>
        &lineCatalogCorrespondingNominalAmp,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti) const {

  for (Int32 k = 0; k != lineCatalogCorrespondingNominalAmp.size(); ++k) {
    Log.LogDebug(Formatter() << "log linesCorrespondingNominalAmp for "
                             << m_lineRatioCatalogs[k].getName());
    for (const auto &[id, line] : lineCatalogCorrespondingNominalAmp[k]) {
      Float64 ebv = enableISMCalzetti
                        ? ismCorrectionCalzetti->GetEbmvValue(GetIsmIndex(k))
                        : NAN;
      Float64 nomAmp = lineCatalogCorrespondingNominalAmp[k].at(id);
      std::string const &lineName = GetCatalog(k).GetList().at(id).GetName();
      Log.LogDebug("    CatalogsTplRatio - "
                   "linesCorrespondingNominalAmp, "
                   "iCatalog=%d, LineId=%d with name=%s, ebv=%f: "
                   "NominalAmpFound = "
                   "%e",
                   k, id, lineName.c_str(), ebv, nomAmp);
    }
  }
}

/**
 * \brief Calculates the best fit between the linemodel fitted amplitudes
 *and the tplRatio catalogs: (for lm-rigidity=tplcorr)
 *
 **/
Float64 CLineCatalogsTplRatio::GetBestFit(const CLineMap &restLineList,
                                          const TFloat64Map &fittedAmplitudes,
                                          const TFloat64Map &fittedErrors,
                                          TFloat64Map &amplitudesCorrected,
                                          std::string &bestTplName) const {

  Float64 coeffMin = -1.0;
  TFloat64Map bestFitAmplitudes;

  auto negativefittedValues = [&fittedAmplitudes,
                               &fittedErrors](Int32 line_id) {
    return (fittedAmplitudes.at(line_id) < 0 || fittedErrors.at(line_id) <= 0);
  };

  // mask is unique and independent from catalogs
  TBoolMap mask;
  for (auto const &[iLine, _] : fittedAmplitudes)
    mask[iLine] = negativefittedValues(iLine) ? false : true;

  // bestFitAmplitudes.resize(restLineList.size());
  for (const auto &catalog : m_lineRatioCatalogs) {

    TFloat64Map ampsCorrected;
    Float64 fit = getFitForOneCatalog(restLineList, fittedAmplitudes,
                                      fittedErrors, catalog, ampsCorrected);

    if (fit > 0.0 && !std::isnan(fit) && (fit < coeffMin || coeffMin == -1.)) {
      coeffMin = fit;
      bestFitAmplitudes = std::move(ampsCorrected);
      bestTplName = catalog.getName();
    }
  }
  // coeff min normalization
  // coeffMin = sqrt(coeffMin);
  // coeffMin /= 1.0;
  amplitudesCorrected.clear();
  if (coeffMin < 0.)
    return coeffMin;

  // fill the corrected amplitudes vector
  Int32 iTplAmps = 0;
  for (auto const &[id, amp] : bestFitAmplitudes) {
    if (mask.at(id))
      amplitudesCorrected.at(id) = amp;
  }

  return coeffMin;
}

Float64 CLineCatalogsTplRatio::getFitForOneCatalog(
    const CLineMap &restLineList, const TFloat64Map &fittedAmplitudes,
    const TFloat64Map &fittedErrors, const CLineRatioCatalog &catalog,
    TFloat64Map &ampsCorrected) const {

  if (fittedAmplitudes.empty())
    return NAN;

  auto const &lineList = catalog.GetList();
  // create the amplitude float vectors
  TFloat64Map tplratioAmplitudes;

  for (auto const &[lineID, line] : restLineList) {

    auto it = lineList.find(lineID);

    if (it == lineList.end())
      tplratioAmplitudes[lineID] = 0.0;
    else
      tplratioAmplitudes[lineID] = it->second.GetNominalAmplitude();
  }

  // ampsCorrected = TFloat64List(linemodelAmplitudes.size());
  Float64 const fit = computeFitValue(fittedAmplitudes, fittedErrors,
                                      tplratioAmplitudes, ampsCorrected);
  return fit;
}

Float64 CLineCatalogsTplRatio::computeFitValue(
    const TFloat64Map &ampsLM, const TFloat64Map &errLM,
    const TFloat64Map &ampsTPL, TFloat64Map &ampsCorrected) const {

  Float64 N = ampsLM.size();

  // // Normalize AmpsLM first
  // Float64 normalizeCoeff =
  //     std::accumulate(ampsLM.cbegin(), ampsLM.cend(), 0.0,
  //                     [](const Float64 previous, const auto &e) {
  //                       return previous + e.second;
  //                     });
  // std::transform(ampsLM.cbegin(), ampsLM.cend(), ampsLM.begin(),
  //                [normalizeCoeff](auto &e) {
  //                  return std::pair(e.first, e.second / normalizeCoeff);
  //                });
  // for (auto &[_, amp] : ampsLM)
  //   amp /= normalizeCoeff;

  Float64 normalizeCoeff = 1.0;

  // estimate fitting amplitude
  Float64 sumLM = 0.0;
  Float64 sumTPL = 0.0;
  Float64 sumCross = 0.0;
  Float64 sumTPL2 = 0.0;

  for (auto const &[id, _] : ampsLM) {
    Float64 err2 = 1.0 / (errLM.at(id) * errLM.at(id));

    sumCross += ampsLM.at(id) * ampsTPL.at(id) * err2;
    sumTPL2 += ampsTPL.at(id) * ampsTPL.at(id) * err2;

    sumLM += ampsLM.at(id) * err2;
    sumTPL += ampsTPL.at(id) * err2;
  }

  if (sumCross == 0 || sumTPL2 == 0)
    return -1.0;

  Float64 ampl = sumCross / sumTPL2;
  Float64 fit = 0.0;
  for (auto const &[id, _] : ampsLM) {
    // compute the chi2
    Float64 const err2 = 1.0 / (errLM.at(id) * errLM.at(id));
    Float64 const diff = ampsLM.at(id) - ampl * ampsTPL.at(id);
    fit += diff * diff * err2;

    // fill the amps_corrected map
    ampsCorrected.at(id) = ampsTPL.at(id) * ampl * normalizeCoeff;
  }

  return fit;
}

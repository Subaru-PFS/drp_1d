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
#include <algorithm> // std::sort
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/line/catalogsTplRatio.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/linemodel/elementlist.h"
#include "RedshiftLibrary/log/log.h"

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;

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

std::vector<std::vector<TFloat64List>>
CLineCatalogsTplRatio::InitLineCorrespondingAmplitudes(
    const std::vector<TLineModelElementParam_ptr> &LineModelElementsParams,
    Int32 enableISMCalzetti,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti) const {

  std::vector<std::vector<TFloat64List>> lineCatalogCorrespondingNominalAmp(
      GetCatalogsCount(),
      std::vector<TFloat64List>(LineModelElementsParams.size()));

  for (Int32 iCatalog = 0; iCatalog != GetCatalogsCount(); ++iCatalog) {
    auto const &catalog = GetCatalog(iCatalog);
    for (Int32 elt_index = 0; elt_index != LineModelElementsParams.size();
         ++elt_index) {
      auto &elt_param_ptr = LineModelElementsParams[elt_index];
      auto &correspondingNominalAmp =
          lineCatalogCorrespondingNominalAmp[iCatalog][elt_index];
      correspondingNominalAmp.reserve(elt_param_ptr->size());
      for (Int32 line_index = 0; line_index != elt_param_ptr->size();
           ++line_index) {
        auto const line_id = elt_param_ptr->m_Lines[line_index].GetID();
        auto const &catalog_line = catalog.GetList().at(line_id);
        Float64 nominalAmp = catalog_line.GetNominalAmplitude();
        Float64 const restLambda = catalog_line.GetPosition();
        if (enableISMCalzetti && elt_param_ptr->m_isEmission == true) {
          Float64 dustCoeff = ismCorrectionCalzetti->GetDustCoeff(
              catalog.getIsmIndex(), restLambda);
          nominalAmp *= dustCoeff;
        }
        correspondingNominalAmp.push_back(nominalAmp);
      }
    }
  }

  // Now log the linesCorrespondingNominalAmp
  logLineNominalAmp(LineModelElementsParams, enableISMCalzetti,
                    lineCatalogCorrespondingNominalAmp, ismCorrectionCalzetti);

  return lineCatalogCorrespondingNominalAmp;
}

// only logging
void CLineCatalogsTplRatio::logLineNominalAmp(
    const std::vector<TLineModelElementParam_ptr> &LineModelElementsParams,
    bool enableISMCalzetti,
    const std::vector<std::vector<TFloat64List>>
        &lineCatalogCorrespondingNominalAmp,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti) const {

  for (Int32 k = 0; k != lineCatalogCorrespondingNominalAmp.size(); ++k) {
    Log.LogDebug(Formatter() << "log linesCorrespondingNominalAmp for "
                             << m_lineRatioCatalogs[k].getName());
    for (Int32 elt_index = 0;
         elt_index != lineCatalogCorrespondingNominalAmp[k].size();
         ++elt_index) {
      for (Int32 line_index = 0;
           line_index !=
           lineCatalogCorrespondingNominalAmp[k][elt_index].size();
           ++line_index) {
        Float64 ebv = enableISMCalzetti
                          ? ismCorrectionCalzetti->GetEbmvValue(GetIsmIndex(k))
                          : NAN;
        Float64 nomAmp =
            lineCatalogCorrespondingNominalAmp[k][elt_index][line_index];
        auto const &line_id =
            LineModelElementsParams[elt_index]->m_Lines[line_index].GetID();
        std::string const &lineName =
            GetCatalog(k).GetList().at(line_id).GetStrID();
        Log.LogDebug("    CatalogsTplRatio - "
                     "linesCorrespondingNominalAmp, "
                     "iCatalog=%d, iElt=%d, iLine=%d with name=%s, ebv=%f: "
                     "NominalAmpFound = "
                     "%e",
                     k, elt_index, line_index, lineName.c_str(), ebv, nomAmp);
      }
    }
  }
}

/**
 * \brief Calculates the best fit between the linemodel fitted amplitudes
 *and the tplRatio catalogs: (for lm-rigidity=tplcorr)
 *
 **/
Float64 CLineCatalogsTplRatio::GetBestFit(const TInt32List &validLinesIndex,
                                          const TFloat64List &fittedAmplitudes,
                                          const TFloat64List &fittedErrors,
                                          TFloat64List &amplitudesCorrected,
                                          Int32 &bestTplRatio) const {

  Float64 coeffMin = -1.0;
  TFloat64List bestFitAmplitudes;

  auto negativefittedValues = [&fittedAmplitudes, &fittedErrors](Int32 idx) {
    return (fittedAmplitudes[idx] < 0 || fittedErrors[idx] <= 0);
  };

  // mask is unique and independent from catalogs
  TBoolList mask;
  for (Int32 idx = 0; idx != validLinesIndex.size(); ++idx)
    mask.push_back(negativefittedValues(idx) ? false : true);

  for (size_t icat = 0; icat < m_lineRatioCatalogs.size(); ++icat) {
    auto const &catalog = m_lineRatioCatalogs[icat];

    TFloat64List ampsCorrected;
    Float64 fit = getFitForOneCatalog(validLinesIndex, fittedAmplitudes,
                                      fittedErrors, catalog, ampsCorrected);

    if (fit > 0.0 && !std::isnan(fit) && (fit < coeffMin || coeffMin == -1.)) {
      coeffMin = fit;
      bestFitAmplitudes = std::move(ampsCorrected);
      bestTplRatio = icat;
    }
  }
  // coeff min normalization
  // coeffMin = sqrt(coeffMin);
  // coeffMin /= 1.0;
  amplitudesCorrected.clear();
  if (coeffMin < 0.)
    return coeffMin;

  // fill the corrected amplitudes vector
  for (Int32 idx = 0; idx != validLinesIndex.size(); ++idx) {
    if (mask[idx])
      amplitudesCorrected[idx] = bestFitAmplitudes[idx];
  }

  return coeffMin;
}

Float64 CLineCatalogsTplRatio::getFitForOneCatalog(
    const TInt32List &validLinesIndex, const TFloat64List &fittedAmplitudes,
    const TFloat64List &fittedErrors, const CLineRatioCatalog &catalog,
    TFloat64List &ampsCorrected) const {

  if (fittedAmplitudes.empty())
    return NAN;

  // create the amplitude float vectors
  // This assumes lineratios and main linecat have all same keys
  TFloat64List tplratioAmplitudes;
  auto const &lineList = catalog.GetList();
  auto it = lineList.cbegin();
  Int32 lastIdx = 0;
  for (auto idx : validLinesIndex) {
    std::advance(it, idx - lastIdx);
    lastIdx = idx;
    tplratioAmplitudes.push_back(it->second.GetNominalAmplitude());
  }

  Float64 const fit = computeFitValue(fittedAmplitudes, fittedErrors,
                                      tplratioAmplitudes, ampsCorrected);
  return fit;
}

Float64 CLineCatalogsTplRatio::computeFitValue(
    const TFloat64List &ampsLM, const TFloat64List &errLM,
    const TFloat64List &ampsTPL, TFloat64List &ampsCorrected) const {

  Float64 N = ampsLM.size();

  // // Normalize AmpsLM first
  // Float64 normalizeCoeff =
  //     std::accumulate(ampsLM.cbegin(), ampsLM.cend(), 0.0,
  //                     [](const Float64 previous, const auto &e) {
  //                       return previous + e;
  //                     });
  // std::transform(ampsLM.cbegin(), ampsLM.cend(), ampsLM.begin(),
  //                [normalizeCoeff](auto &e) {
  //                  return e / normalizeCoeff);
  //                });
  // for (auto &amp : ampsLM)
  //   amp /= normalizeCoeff;

  Float64 normalizeCoeff = 1.0;

  // estimate fitting amplitude
  Float64 sumLM = 0.0;
  Float64 sumTPL = 0.0;
  Float64 sumCross = 0.0;
  Float64 sumTPL2 = 0.0;

  for (Int32 i = 0; i != ampsLM.size(); ++i) {
    Float64 err2 = 1.0 / (errLM[i] * errLM[i]);

    sumCross += ampsLM[i] * ampsTPL[i] * err2;
    sumTPL2 += ampsTPL[i] * ampsTPL[i] * err2;

    sumLM += ampsLM[i] * err2;
    sumTPL += ampsTPL[i] * err2;
  }

  if (sumCross == 0 || sumTPL2 == 0)
    return -1.0;

  Float64 ampl = sumCross / sumTPL2;
  Float64 fit = 0.0;
  ampsCorrected.assign(ampsLM.size(), NAN);
  for (Int32 i = 0; i != ampsLM.size(); ++i) {
    // compute the chi2
    Float64 const err2 = 1.0 / (errLM[i] * errLM[i]);
    Float64 const diff = ampsLM[i] - ampl * ampsTPL[i];
    fit += diff * diff * err2;

    // fill the amps_corrected map
    ampsCorrected[i] = ampsTPL[i] * ampl * normalizeCoeff;
  }

  return fit;
}

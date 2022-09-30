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
CLineCatalog::TLineVector
CLineCatalogsTplRatio::GetRestLinesList(Int32 index) const {
  Int32 typeFilter = -1;
  Int32 forceFilter = -1;

  CLineCatalog::TLineVector restLineList =
      m_lineRatioCatalogs[index].GetFilteredList(typeFilter, forceFilter);
  return restLineList;
}

Int32 CLineCatalogsTplRatio::GetCatalogsCount() const {
  return m_lineRatioCatalogs.size();
}

TFloat64List CLineCatalogsTplRatio::getCatalogsPriors() const {
  TFloat64List catalogsPriors;
  if (catalogsPriors.empty()) {
    for (CLineRatioCatalog cat : m_lineRatioCatalogs)
      catalogsPriors.push_back(cat.getPrior());
  }
  return catalogsPriors;
}

std::string CLineCatalogsTplRatio::GetCatalogName(Int32 idx) const {
  return m_lineRatioCatalogs[idx].getName();
}

Int32 CLineCatalogsTplRatio::GetIsmIndex(Int32 idx) const {
  return m_lineRatioCatalogs[idx].getIsmIndex();
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
    const CLineModelElementList &LineModelElementList, Int32 enableISMCalzetti,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti,
    Float64 nsigmasupport) const {

  std::vector<std::vector<TFloat64List>>
      lineCatalogLinesCorrespondingNominalAmp;

  lineCatalogLinesCorrespondingNominalAmp.reserve(LineModelElementList.size());
  const Int32 catalogCount = GetCatalogsCount();
  for (const auto &elt : LineModelElementList) {
    // first set all corresponding amplitudes to 0.0;
    lineCatalogLinesCorrespondingNominalAmp.push_back(std::vector<TFloat64List>(
        catalogCount, TFloat64List(elt->GetSize(), 0.0)));

    // now set the non-zero amp correspondences
    for (Int32 iCatalog = 0; iCatalog < catalogCount; iCatalog++) {
      const CLineCatalog::TLineVector lineList =
          m_lineRatioCatalogs[iCatalog].GetList();

      for (const auto &currentline : lineList) {
        const std::string &currentLineName = currentline.GetName();
        Float64 nominalAmp = currentline.GetNominalAmplitude();
        Float64 restLambda = currentline.GetPosition();
        if (enableISMCalzetti) {
          Float64 dustCoeff = ismCorrectionCalzetti->GetDustCoeff(
              GetIsmIndex(iCatalog), restLambda);
          nominalAmp *= dustCoeff;
        }
        // find line in the elementList and fill with nominalAmp
        auto it = std::find_if(elt->m_Lines.cbegin(), elt->m_Lines.cend(),
                               [&currentLineName](const CLine &line) {
                                 return line.GetName() == currentLineName;
                               });
        if (it != elt->m_Lines.end()) {
          Int32 idx = it - elt->m_Lines.begin();
          lineCatalogLinesCorrespondingNominalAmp.back()[iCatalog][idx] =
              nominalAmp;
        }
      }
    }
  }
  // Now log the linesCorrespondingNominalAmp
  logLineNominalAmp(LineModelElementList, enableISMCalzetti,
                    lineCatalogLinesCorrespondingNominalAmp,
                    ismCorrectionCalzetti);

  return lineCatalogLinesCorrespondingNominalAmp;
}

// only logging
void CLineCatalogsTplRatio::logLineNominalAmp(
    const CLineModelElementList &LineModelElementList, bool enableISMCalzetti,
    const std::vector<std::vector<TFloat64List>>
        &lineCatalogLinesCorrespondingNominalAmp,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        &ismCorrectionCalzetti) const {

  for (Int32 iElt = 0; iElt < LineModelElementList.size(); iElt++) {
    for (Int32 k = 0; k < GetCatalogsCount(); k++) {
      /*Log.LogDebug(Formatter() << "log linesCorrespondingNominalAmp for "
                               << m_lineRatioCatalogs[k].getName());*/
      Int32 nLines = LineModelElementList[iElt]->GetSize();
      for (Int32 j = 0; j < nLines; j++) {
        Float64 ebv = enableISMCalzetti
                          ? ismCorrectionCalzetti->GetEbmvValue(GetIsmIndex(k))
                          : NAN;
        Float64 nomAmp = lineCatalogLinesCorrespondingNominalAmp[iElt][k][j];
        std::string lineName = LineModelElementList[iElt]->m_Lines[j].GetName();
        /*Log.LogDebug("    CatalogsTplRatio - "
                     "linesCorrespondingNominalAmp iElt=%d, "
                     "iCatalog=%d, iLine=%d with name=%s, ebv=%f: "
                     "NominalAmpFound = "
                     "%e",
                     iElt, k, j, lineName.c_str(), ebv, nomAmp);*/
      }
    }
  }
}

const CLineCatalog &CLineCatalogsTplRatio::GetCatalog(Int32 iCatalog) const {
  return m_lineRatioCatalogs[iCatalog];
}

/**
 * \brief Calculates the best fit between the linemodel fitted amplitudes
 *and the tplRatio catalogs: (for lm-rigidity=tplcorr)
 *
 **/
Float64 CLineCatalogsTplRatio::GetBestFit(
    const CLineCatalog::TLineVector &restLineList,
    const TFloat64List &fittedAmplitudes, const TFloat64List &fittedErrors,
    TFloat64List &amplitudesCorrected, std::string &bestTplName) const {

  Float64 coeffMin = -1.0;
  TFloat64List bestFitAmplitudes;

  auto negativefittedValues = [&fittedAmplitudes, &fittedErrors](Int32 idx) {
    return (fittedAmplitudes[idx] < 0 || fittedErrors[idx] <= 0);
  };

  // mask is unique and independent from catalogs
  TInt32List mask(fittedAmplitudes.size());
  for (Int32 iRestLine = 0; iRestLine < restLineList.size(); iRestLine++) {
    if (negativefittedValues(iRestLine)) {
      mask[iRestLine] = 0;
      continue;
    }
    mask[iRestLine] = 1;
  }

  // bestFitAmplitudes.resize(restLineList.size());
  for (const auto &catalog : m_lineRatioCatalogs) {

    TFloat64List ampsCorrected;
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
  amplitudesCorrected =
      TFloat64List(restLineList.size(), -1.); // NAN plutot que -1?
  if (coeffMin < 0.)                          // no need to check mask value
    return coeffMin;

  // fill the corrected amplitudes vector
  Int32 iTplAmps = 0;
  for (Int32 iRestLine = 0; iRestLine < restLineList.size(); iRestLine++) {
    if (mask[iRestLine] == 1) {
      amplitudesCorrected[iRestLine] = bestFitAmplitudes[iTplAmps];
      iTplAmps++;
    }
  }

  return coeffMin;
}

Float64 CLineCatalogsTplRatio::getFitForOneCatalog(
    const CLineCatalog::TLineVector &restLineList,
    const TFloat64List &fittedAmplitudes, const TFloat64List &fittedErrors,
    const CLineRatioCatalog &catalog, TFloat64List &ampsCorrected) const {

  CLineCatalog::TLineVector lineList = catalog.GetList();
  // create the amplitude float vectors
  TFloat64List tplratioAmplitudes;
  TFloat64List linemodelAmplitudes;
  TFloat64List linemodelErrors;
  for (Int32 iRestLine = 0; iRestLine < restLineList.size(); iRestLine++) {
    std::string restLineName = restLineList[iRestLine].GetName();

    auto it = std::find_if(lineList.cbegin(), lineList.cend(),
                           [restLineName](const CLine &line) {
                             return line.GetName() == restLineName;
                           });
    if (it == lineList.end())
      tplratioAmplitudes.push_back(0.0);
    else {
      Int32 tplRatioIdx = it - lineList.begin();
      Float64 amp = lineList[tplRatioIdx].GetNominalAmplitude();
      tplratioAmplitudes.push_back(amp);
    }
    linemodelAmplitudes.push_back(fittedAmplitudes[iRestLine]);
    linemodelErrors.push_back(fittedErrors[iRestLine]);
  }

  if (!linemodelAmplitudes.size())
    return NAN;

  ampsCorrected = TFloat64List(linemodelAmplitudes.size());
  Float64 fit = computeFitValue(linemodelAmplitudes, linemodelErrors,
                                tplratioAmplitudes, ampsCorrected);
  return fit;
}

Float64 CLineCatalogsTplRatio::computeFitValue(
    const TFloat64List &ampsLM, const TFloat64List &errLM,
    const TFloat64List &ampsTPL, TFloat64List &ampsCorrected) const {

  Float64 N = ampsLM.size();
  //    // Normalize AmpsLM first
  //    Float64 normalizeCoeff = 0.0;
  //    for(Int32 k=0; k<N; k++)
  //    {
  //        normalizeCoeff+= ampsLM[k];
  //    }
  //    for(Int32 k=0; k<N; k++)
  //    {
  //        ampsLM[k] = ampsLM[k]/normalizeCoeff;
  //    }
  Float64 normalizeCoeff = 1.0;

  // estimate fitting amplitude
  Float64 sumLM = 0.0;
  Float64 sumTPL = 0.0;
  Float64 sumCross = 0.0;
  Float64 sumTPL2 = 0.0;

  for (Int32 k = 0; k < N; k++) {
    Float64 err2 = 1.0 / (errLM[k] * errLM[k]);

    sumCross += ampsLM[k] * ampsTPL[k] * err2;
    sumTPL2 += ampsTPL[k] * ampsTPL[k] * err2;

    sumLM += ampsLM[k] * err2;
    sumTPL += ampsTPL[k] * err2;
  }
  //    if ( sumLM==0 || sumTPL==0 )
  //    {
  //        return -1.0;
  //    }
  //    Float64 ampl = sumLM / sumTPL;
  if (sumCross == 0 || sumTPL2 == 0) {
    return -1.0;
  }
  Float64 ampl = sumCross / sumTPL2;

  Float64 fit = 0.0;
  Float64 diff;
  for (Int32 k = 0; k < N; k++) {
    Float64 err2 = 1.0 / (errLM[k] * errLM[k]);
    diff = ampsLM[k] - ampl * ampsTPL[k];
    fit += diff * diff * err2;
  }

  // fill the amps_corrected vector
  for (Int32 k = 0; k < N; k++) {
    ampsCorrected[k] = ampsTPL[k] * ampl * normalizeCoeff;
  }

  return fit;
}

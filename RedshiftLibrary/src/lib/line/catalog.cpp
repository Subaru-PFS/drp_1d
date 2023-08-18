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
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/line/linetags.h"
#include "RedshiftLibrary/log/log.h"

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <algorithm> // std::sort
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace NSEpic;
using namespace std;
using namespace boost;
using namespace boost::filesystem;

template <typename TLine>
std::vector<TLine>
CLineCatalogBase<TLine>::GetFilteredList(CLine::EType typeFilter,
                                         CLine::EForce forceFilter) const {

  TLineVector filteredList;

  for (Int32 i = 0; i < m_List.size(); i++) {
    if (typeFilter == CLine::EType::nType_All ||
        typeFilter == m_List[i].GetType()) {
      if (forceFilter == CLine::EForce::nForce_All ||
          forceFilter == m_List[i].GetForce()) {
        filteredList.push_back(m_List[i]);
      }
    }
  }

  /*   for (const auto &it : m_List) {
      const auto &str_Id = it.first;
      const auto &line = it.second;
      if ((typeFilter == CLine::EType::nType_All ||
           typeFilter == line.GetType()) &&
          (forceFilter == CLine::EForce::nForce_All ||
           forceFilter == line.GetForce())) {
        filteredList[str_Id] = line;
      }
    }
   */

  return filteredList;
}

template <typename TLine>
std::vector<TLine>
CLineCatalogBase<TLine>::GetFilteredList(const std::string &typeFilter,
                                         const std::string &forceFilter) const {
  const auto &etypeFilter = CLine::string2Type(typeFilter);
  const auto &eforceFilter = CLine::string2Force(forceFilter);

  return GetFilteredList(etypeFilter, eforceFilter);
}

template <typename TLine>
const std::vector<std::vector<TLine>>
CLineCatalogBase<TLine>::ConvertToGroupList(
    const std::vector<TLine> &filteredList) {

  std::vector<TLineVector> fullList;

  TStringList tags;
  for (int i = 0; i < filteredList.size(); i++) {
    if (filteredList[i].GetGroupName() != undefStr) {
      tags.push_back(filteredList[i].GetGroupName());
    }
  }

  // create the group tag set by removing duplicates
  std::sort(tags.begin(), tags.end());
  tags.erase(std::unique(tags.begin(), tags.end()), tags.end());

  // get all group tags
  for (Int32 itag = 0; itag < tags.size(); itag++) {
    std::vector<TLine> taggedGroupList;
    for (int i = 0; i < filteredList.size(); i++) {
      std::string group = filteredList[i].GetGroupName();
      if (group == tags[itag]) {
        taggedGroupList.push_back(filteredList[i]);
      }
    }
    fullList.push_back(taggedGroupList);
  }
  // add the non grouped lines
  for (int i = 0; i < filteredList.size(); i++) {
    std::string group = filteredList[i].GetGroupName();
    if (group == undefStr) {
      std::vector<TLine> taggedGroupList;
      taggedGroupList.push_back(filteredList[i]);
      fullList.push_back(taggedGroupList);
    }
  }

  /*
    for (const auto &it : filteredList) {
      const auto &line = it.second;
      const auto &str_Id = it.first;
      const auto &group_name = line.GetGroupName();
      if (group_name == undefStr)
        // non grouped lines are added in dedicated maps (one element)
        filteredList[str_Id] = std::vector<TLine>{{str_Id, line}};
      else
        filteredList[group_name][str_Id] = line;
    }
  */
  return fullList;
}

template <typename TLine> void CLineCatalogBase<TLine>::Add(const TLine &r) {

  typename TLineVector::iterator it;
  for (it = m_List.begin(); it != m_List.end(); ++it) {
    if ((*it).GetPosition() == r.GetPosition() &&
        (*it).GetType() == r.GetType())
      THROWG(INTERNAL_ERROR,
             Formatter() << "line with name " << r.GetName()
                         << " already exists, position=" << r.GetPosition()
                         << " type=" << r.GetTypeString());
  }

  m_List.push_back(r);

  /*   const auto &str_Id = r.GetStrID();
    // verify unicity
    if (m_List.find(str_Id) != m_List.end())
      THROWG(INTERNAL_ERROR, Formatter() << "Duplicated Ids: line " << str_Id
                                          << "is already in the catalog");

    // insert new line
    m_List[str_Id] = r;
   */
}

void CLineCatalog::AddLineFromParams(
    const std::string &name, const Float64 &position, const std::string &type,
    const std::string &force, const std::string &profileName,
    const TAsymParams &asymParams, const std::string &groupName,
    const Float64 &nominalAmplitude, const std::string &velocityGroup,
    const Float64 &velocityOffset, const bool &enableVelocityFit,
    const Int32 &id, const std::string &str_id,
    const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrection) {

  const auto &etype = CLine::string2Type(type);
  const auto &eforce = CLine::string2Force(force);

  TAsymParams _asymParams = {1., 4.5, 0.};
  TAsymParams _asymFitParams = {2., 2., 0.};
  std::unique_ptr<CLineProfile> profile;

  if (profileName.find("ASYMFIXED") != std::string::npos)
    profile = std::unique_ptr<CLineProfileASYM>(
        new CLineProfileASYM(m_nSigmaSupport, asymParams, "mean"));
  else if (profileName == "SYM")
    profile =
        std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM(m_nSigmaSupport));
  else if (profileName == "LOR")
    profile =
        std::unique_ptr<CLineProfileLOR>(new CLineProfileLOR(m_nSigmaSupport));
  else if (profileName == "ASYM")
    profile = std::unique_ptr<CLineProfileASYM>(
        new CLineProfileASYM(m_nSigmaSupport, _asymParams, "none"));
  else if (profileName == "ASYMFIT")
    profile = std::unique_ptr<CLineProfileASYMFIT>(
        new CLineProfileASYMFIT(m_nSigmaSupport, _asymFitParams, "mean"));
  else if (profileName == "SYMIGM")
    profile = std::unique_ptr<CLineProfileSYMIGM>(
        new CLineProfileSYMIGM(igmcorrection, m_nSigmaSupport));
  else {
    THROWG(INTERNAL_ERROR, Formatter() << "Profile name " << profileName
                                       << " is no recognized.");
  }

  Add(CLine(name, position, etype, std::move(profile), eforce, velocityOffset,
            enableVelocityFit, groupName, nominalAmplitude, velocityGroup, id,
            str_id));
}

void CLineDetectedCatalog::Sort() { sort(m_List.begin(), m_List.end()); }

template <typename TLine>
void CLineCatalogBase<TLine>::setLineAmplitude(
    const std::string &str_id, const Float64 &nominalAmplitude) {
  typename TLineVector::iterator it;
  for (it = m_List.begin(); it != m_List.end(); ++it) {
    if (it->GetStrID() == str_id)
      return it->setNominalAmplitude(nominalAmplitude);
  }
  THROWG(INTERNAL_ERROR, Formatter() << " Line with id " << str_id
                                     << " does not exist in catalog");

  /*   const auto &search = m_List.find(str_id);
    if (search != m_List.end())
      search->second.setNominalAmplitude(nominalAmplitude);
    else
      THROWG(INTERNAL_ERROR, Formatter() << " Line with id " << str_id
                                         << " does not exist in catalog");
   */
}

/**
 * @brief only called for tplRatio catalogs
 * hereinafter we consider that only one ASYMFIT profile exists in linecatalog
 * and that it corresponds to Lya
 * @param profile
 * @param params
 */
template <typename TLine>
void CLineCatalogBase<TLine>::setAsymProfileAndParams(
    const std::string &profile, TAsymParams params) {

  typename TLineVector::iterator it;
  for (it = m_List.begin(); it != m_List.end(); ++it) {
    if (it->GetProfile()->isAsymFit() || it->GetProfile()->isAsymFixed())
      return it->setProfileAndParams(profile, params, m_nSigmaSupport);
  }

  /*   for (const auto it : m_List) {
      const auto &line = it.second;
      if (line.GetProfile()->isAsymFit() || line.GetProfile()->isAsymFixed())
        line.setProfileAndParams(profile, params, m_nSigmaSupport);
    }
   */
}

template <typename TLine>
void CLineCatalogBase<TLine>::convertLineProfiles2SYMIGM(
    const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> &igmcorrection) {

  linetags ltags;
  typename TLineVector::iterator it;
  for (it = m_List.begin(); it != m_List.end(); ++it) {
    if (!it->IsEmission())
      continue;
    if (it->GetName() == ltags.lya_em || it->GetPosition() < RESTLAMBDA_LYA)
      it->setProfileAndParams("SYMIGM", TAsymParams(), m_nSigmaSupport,
                              igmcorrection);
  }

  /*   for (const auto &it : m_List) {
      const auto &line = it.second;
      if (!line.IsEmission())
        continue;
      if (line.GetName() == linetags::lya_em ||
          line.GetPosition() < RESTLAMBDA_LYA)
        line.setProfileAndParams("SYMIGM", TAsymParams(), m_nSigmaSupport,
                                 igmcorrection);
    } */
}

template class CLineCatalogBase<CLine>;
template class CLineCatalogBase<CLineDetected>;

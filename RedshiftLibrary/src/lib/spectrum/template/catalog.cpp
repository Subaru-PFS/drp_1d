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
#include <string>

#include <boost/filesystem.hpp>

#include "RedshiftLibrary/continuum/irregularsamplingmedian.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/spectrum/template/template.h"

using namespace NSEpic;
using namespace std;
using namespace boost::filesystem;

/**
 * Variable instantiator constructor.
 */
CTemplateCatalog::CTemplateCatalog(bool sampling) { m_logsampling = sampling; }

TTemplateConstRefList
CTemplateCatalog::const_TTemplateRefList_cast(const TTemplateRefList &list) {
  TTemplateConstRefList const_list;
  for (auto tpl : list)
    const_list.push_back(tpl);

  return const_list;
}

/**
 * Returns a list containing all templates as enumerated in the categoryList
 * input.
 */
TTemplateRefList
CTemplateCatalog::GetTemplateList_(const TStringList &categoryList) const {

  return GetTemplateList_(categoryList, m_orthogonal, m_logsampling);
}

/**
 * Returns a list containing all templates as enumerated in the categoryList
 * input.
 */
TTemplateRefList
CTemplateCatalog::GetTemplateList_(const TStringList &categoryList,
                                   bool opt_ortho, bool opt_logsampling) const {
  TTemplateRefList list;

  for (const auto &cat : categoryList) {
    const auto &map = GetList(opt_ortho, opt_logsampling);
    if (!map.count(cat))
      continue;
    for (const auto &tpl : map.at(cat))
      list.push_back(tpl);
  }
  return list;
}

std::shared_ptr<const CTemplate>
CTemplateCatalog::GetTemplateByName(const TStringList &tplCategoryList,
                                    const std::string tplName, bool opt_ortho,
                                    bool opt_logsampling) const {
  for (const auto &cat : tplCategoryList) {
    const auto &map = GetList(opt_ortho, opt_logsampling);
    if (!map.count(cat))
      continue;
    for (const auto &tpl : map.at(cat))
      if (tpl->GetName() == tplName)
        return tpl;
  }
  THROWG(ErrorCode::INTERNAL_ERROR,
         Formatter() << "Could not find template with name " << tplName);
}

std::shared_ptr<const CTemplate>
CTemplateCatalog::GetTemplateByName(const TStringList &tplCategoryList,
                                    const std::string tplName) const {
  return GetTemplateByName(tplCategoryList, tplName, m_orthogonal,
                           m_logsampling);
}

/**
 * Get a list of strings with the contents of m_List.
 */
TStringList CTemplateCatalog::GetCategoryList() const {
  TStringList l;
  for (const auto &it : GetList()) {
    l.push_back(it.first);
  }
  return l;
}

/**
 * Returns the size of the category entry in m_List.
 */
Int32 CTemplateCatalog::GetTemplateCount(const std::string &category) const {
  return GetTemplateCount(category, m_orthogonal, m_logsampling);
}

Int32 CTemplateCatalog::GetTemplateCount(const std::string &category,
                                         bool opt_ortho,
                                         bool opt_logsampling) const {
  if (!GetList(opt_ortho, opt_logsampling).count(category))
    return 0;
  return GetList(opt_ortho, opt_logsampling).at(category).size();
}

Int32 CTemplateCatalog::GetNonNullTemplateCount(const std::string &category,
                                                bool opt_ortho,
                                                bool opt_logsampling) const {
  if (!GetTemplateCount(category, opt_ortho, opt_logsampling))
    return 0;

  const TTemplatesRefDict &tplList = GetList(opt_ortho, opt_logsampling);
  Int32 count = 0;
  for (auto tpl : tplList.at(category))
    if (tpl)
      count++;

  return count;
}

/**
 * Adds the input to the list of templates, under its category. If the input
 * doesn't have a category, function returns false. Also computes the template
 * without continuum and adds it to the list of templates without continuum.
 * Returns true.
 * @sampling here is relevant here especially for when rebinning
 * otherwise we will have to SetSampling ("lin") when we want to read the
 * non-log template, and then call SetSampling("log") when we want to add the
 * log tpl
 */
void CTemplateCatalog::Add(const std::shared_ptr<CTemplate> &r) {
  if (r->GetCategory().empty())
    THROWG(ErrorCode::INTERNAL_ERROR, "Template has no category");

  GetList()[r->GetCategory()].push_back(r);
}

void CTemplateCatalog::SetTemplate(const std::shared_ptr<CTemplate> &tpl,
                                   Int32 i) {
  const std::string category = tpl->GetCategory();

  // first clear given template
  ClearTemplates(category, m_orthogonal, m_logsampling, i);

  // assign new template
  if (!GetTemplateCount(category))
    Add(tpl);
  else
    GetList().at(category).at(i) = tpl;
}

// clear all templates in a category
void CTemplateCatalog::ClearTemplateList(const std::string &category,
                                         bool opt_ortho, bool opt_logsampling) {
  ClearTemplates(category, opt_ortho, opt_logsampling, 0, true);
}

// clear one template in a category
void CTemplateCatalog::ClearTemplates(const std::string &category,
                                      bool opt_ortho, bool opt_logsampling,
                                      Int32 i, bool alltemplates) {
  if (!GetTemplateCount(category, opt_ortho, opt_logsampling))
    return;
  TTemplatesRefDict &tplList = GetList(opt_ortho, opt_logsampling);

  if (alltemplates) {
    tplList.at(category).clear();
    if (opt_ortho) { // reset LSFWidth if the 2 ortho list are empty
      auto &otherOrthoList = GetList(opt_ortho, !opt_logsampling);
      if (otherOrthoList.find(category) == otherOrthoList.end() ||
          otherOrthoList.at(category).empty())
        m_ortho_LSFWidth = NAN;
    }
  } else {
    tplList.at(category).at(i) = nullptr; // clear one element in the list

    // clear the list if all nullptr
    if (!GetNonNullTemplateCount(category, opt_ortho, opt_logsampling))
      ClearTemplates(category, opt_ortho, opt_logsampling, 0, true);
  }

  // Other templates clearing rules
  if (!opt_ortho) // if not ortho clear associated ortho
    ClearTemplates(category, 1, opt_logsampling, i, alltemplates);
  if (!opt_logsampling && !opt_ortho) { // if original, clear log-sampled
    ClearTemplates(category, 0, 1, i,
                   alltemplates); // clear log-sampled (will clear log-sampled
                                  // ortho recursively)
  }
}

void CTemplateCatalog::InitContinuumRemoval(
    const std::shared_ptr<const CParameterStore> &parameterStore) {
  std::string ContinuumRemovalMethod = parameterStore->Get<std::string>(
      "templateCatalog.continuumRemoval.method");
  Float64 medianKernelWidth = parameterStore->Get<Float64>(
      "templateCatalog.continuumRemoval.medianKernelWidth");
  bool medianEvenReflection = parameterStore->Get<bool>(
      "templateCatalog.continuumRemoval.medianEvenReflection");
  for (auto it : GetList(0, 0)) {
    const TTemplateRefList &TplList = it.second;
    for (auto tpl : TplList) {
      tpl->SetContinuumEstimationMethod(ContinuumRemovalMethod);
      tpl->SetMedianWinsize(medianKernelWidth);
      tpl->SetMedianEvenReflection(medianEvenReflection);
    }
  }
}

void CTemplateCatalog::SetIsmIgmCorrection(
    const std::shared_ptr<const CParameterStore> &parameterStore,
    const std::shared_ptr<const CSpectrumFluxCorrectionMeiksin>
        igmCorrectionMeiksin,
    const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>
        ismCorrectionCalzetti) {

  for (const auto &it : GetList(0, 0)) {
    const auto &cat = it.first;
    if (parameterStore->HasTplIsmExtinction(cat) && !it.second.empty() &&
        it.second[0]->CalzettiInitFailed()) {
      for (const auto &tpl : it.second) {
        tpl->m_ismCorrectionCalzetti = ismCorrectionCalzetti;
      }
    }
    if (parameterStore->HasTplIgmExtinction(cat) && !it.second.empty() &&
        it.second[0]->MeiksinInitFailed()) {
      for (const auto &tpl : it.second) {
        tpl->m_igmCorrectionMeiksin = igmCorrectionMeiksin;
      }
    }
  }
}

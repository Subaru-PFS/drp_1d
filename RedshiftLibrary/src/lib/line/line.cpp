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

#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;

const std::map<CLine::EType, std::string> CLine::ETypeString = {
    {EType::nType_Absorption, "Abs"},
    {EType::nType_Emission, "Em"},
    {EType::nType_All, "All"}};

const std::map<CLine::EForce, std::string> CLine::EForceString = {
    {EForce::nForce_Weak, "Weak"}, {EForce::nForce_Strong, "Strong"}};

CLine::EType CLine::string2Type(std::string const &stype) {
  if (stype == "no")
    return EType::nType_All;
  auto const ctype = stype.front();
  if (ctype == 'E')
    return EType::nType_Emission;
  if (ctype == 'A')
    return EType::nType_Absorption;

  THROWG(BAD_LINE_TYPE, Formatter()
                            << "Bad line type, should be in {A,E} : " << ctype);
}

CLine::EForce CLine::string2Force(std::string const &sforce) {
  if (sforce == "no")
    return EForce::nForce_All;
  auto const cforce = sforce.front();
  if (cforce == 'W')
    return EForce::nForce_Weak;
  if (cforce == 'S')
    return EForce::nForce_Strong;

  THROWG(BAD_LINE_FORCE,
         Formatter() << "Bad line force, should be in {S,W} : " << cforce);
}

CLine::CLine(const string &name, Float64 pos, EType type,
             CLineProfile_ptr &&profile, EForce force, Float64 velocityOffset,
             bool enableVelocityOffsetFitting, const std::string &groupName,
             Float64 nominalAmp, const string &velGroupName, Int32 id,
             const std::string &str_id)
    : m_Name(name), m_Pos(pos), m_Type(type), m_Force(force),
      m_Profile(std::move(profile)), m_GroupName(groupName),
      m_NominalAmplitude(nominalAmp), m_VelGroupName(velGroupName), m_id(id),
      m_Offset(velocityOffset), m_OffsetFit(enableVelocityOffsetFitting),
      m_strID(str_id) {}

CLine::CLine(const CLine &other)
    : m_Type(other.m_Type),
      m_Profile(other.m_Profile->Clone()), // deep copy for m_Profile
      m_Force(other.m_Force), m_Pos(other.m_Pos), m_Offset(other.m_Offset),
      m_Name(other.m_Name), m_GroupName(other.m_GroupName),
      m_NominalAmplitude(other.m_NominalAmplitude),
      m_OffsetFit(other.m_OffsetFit), m_VelGroupName(other.m_VelGroupName),
      m_id(other.m_id), m_strID(other.m_strID) {}

CLine &CLine::operator=(const CLine &other) {
  m_id = other.m_id;
  m_Type = other.m_Type;
  m_Profile = other.m_Profile->Clone(); // deep copy for m_Profile
  m_Force = other.m_Force;
  m_Pos = other.m_Pos;
  m_Offset = other.m_Offset;

  m_Name = other.m_Name;

  m_GroupName = other.m_GroupName;
  m_NominalAmplitude = other.m_NominalAmplitude;

  m_OffsetFit = other.m_OffsetFit;

  m_VelGroupName = other.m_VelGroupName;

  m_strID = other.m_strID;

  return *this;
}

void CLine::SetAsymParams(const TAsymParams &asymParams) {
  if (!m_Profile)
    THROWG(INTERNAL_ERROR, "lineprofile is not initialized");
  m_Profile->SetAsymParams(asymParams);
}

void CLine::SetSymIgmParams(const TSymIgmParams &params) {
  if (!m_Profile)
    THROWG(INTERNAL_ERROR, "lineprofile is not initialized");
  m_Profile->SetSymIgmParams(params);
}

void CLine::SetSymIgmFit(bool val) {
  if (!m_Profile)
    THROWG(INTERNAL_ERROR, "lineprofile is not initialized");
  m_Profile->SetSymIgmFit(val);
}

void CLine::setProfileAndParams(
    const std::string &profileName, const TAsymParams &asymParams,
    Float64 nSigmaSupport,
    const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>
        &igmcorrectionMeiksin) {
  if (profileName.find("ASYMFIXED") != std::string::npos)
    m_Profile.reset(new CLineProfileASYM(nSigmaSupport, asymParams, "mean"));
  else if (profileName == "SYM")
    m_Profile.reset(new CLineProfileSYM(nSigmaSupport));
  else if (profileName == "LOR")
    m_Profile.reset(new CLineProfileLOR(nSigmaSupport));
  else if (profileName == "ASYM")
    m_Profile.reset(new CLineProfileASYM(nSigmaSupport, asymParams, "none"));
  else if (profileName == "ASYMFIT")
    m_Profile.reset(new CLineProfileASYMFIT(nSigmaSupport, asymParams, "mean"));
  else if (profileName == "SYMIGM")
    m_Profile.reset(
        new CLineProfileSYMIGM(igmcorrectionMeiksin, nSigmaSupport));
  else
    THROWG(INTERNAL_ERROR, Formatter() << "Profile name " << profileName
                                       << " is no recognized.");
}

void CLine::resetAsymFitParams() {
  if (!m_Profile)
    THROWG(INTERNAL_ERROR, "lineprofile is not initialized");
  m_Profile->resetParams();
}

TAsymParams CLine::GetAsymParams() const {
  if (!m_Profile)
    THROWG(INTERNAL_ERROR, "lineprofile is not initialized");
  return m_Profile->GetAsymParams();
}

TSymIgmParams CLine::GetSymIgmParams() const {
  if (!m_Profile)
    THROWG(INTERNAL_ERROR, "lineprofile is not initialized");
  return m_Profile->GetSymIgmParams();
}

const CLineProfile_ptr &CLine::GetProfile() const {
  if (!m_Profile)
    THROWG(INTERNAL_ERROR, "Current line does not have a profile");
  return m_Profile;
}
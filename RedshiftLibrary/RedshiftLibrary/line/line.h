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
#ifndef _REDSHIFT_LINE_LINE_
#define _REDSHIFT_LINE_LINE_

#include <string>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/line/lineprofile.h"
#include "RedshiftLibrary/line/lineprofileASYM.h"
#include "RedshiftLibrary/line/lineprofileASYMFIT.h"
#include "RedshiftLibrary/line/lineprofileLOR.h"
#include "RedshiftLibrary/line/lineprofileSYM.h"
#include "RedshiftLibrary/line/lineprofileSYMIGM.h"

namespace NSEpic {

/**
 * \ingroup Redshift
 * Represent a Single Line.
 */
class CSpectrumFluxCorrectionMeiksin;
class CLine {

public:
  enum class EType {
    nType_Absorption = 1,
    nType_Emission = 2,
    nType_All = 3,
  };

  static const std::map<EType, std::string> ETypeString;

  enum class EForce {
    nForce_Weak = 1,
    nForce_Strong = 2,
    nForce_All = 3,
  };

  static const std::map<EForce, std::string> EForceString;

  CLine() = default;
  CLine(const std::string &name, Float64 pos, EType type,
        CLineProfile_ptr &&profile, EForce force, Float64 velocityOffset,
        bool enableVelocityOffsetFitting, const std::string &groupName,
        Float64 nominalAmp, const std::string &velGroupName, Int32 id,
        const std::string &str_id);
  virtual ~CLine() = default;
  CLine(const CLine &other);
  CLine(CLine &&other) = default;
  CLine &operator=(const CLine &other);
  CLine &operator=(CLine &&other) = default;

  bool operator!=(const CLine &rhs) const { return (m_id != rhs.m_id); };

  Int32 GetID() const { return m_id; };
  bool IsStrong() const { return m_Force == EForce::nForce_Strong; }
  bool IsWeak() const { return m_Force == EForce::nForce_Weak; }
  bool IsEmission() const { return m_Type == EType::nType_Emission; };
  EForce GetForce() const { return m_Force; };
  EType GetType() const { return m_Type; };
  std::string GetTypeString() const { return ETypeString.at(m_Type); };
  std::string GetForceString() const { return EForceString.at(m_Force); };
  const CLineProfile_ptr &GetProfile() const;
  void SetProfile(CLineProfile_ptr &&profile) {
    m_Profile = std::move(profile);
  };
  Float64 GetPosition() const { return m_Pos; };
  Float64 GetOffset() const { return m_Offset; };
  bool IsOffsetFitEnabled() const { return m_OffsetFit; };
  bool EnableOffsetFit() {
    m_OffsetFit = true;
    return true;
  };
  bool DisableOffsetFit() {
    m_OffsetFit = false;
    return true;
  };

  TAsymParams GetAsymParams() const;
  TSymIgmParams GetSymIgmParams() const;
  void SetAsymParams(const TAsymParams &asymParams);
  void SetSymIgmParams(const TSymIgmParams &params);
  void SetSymIgmFit(bool val = true);
  void setProfileAndParams(const std::string &profileName,
                           const TAsymParams &asymParams, Float64 nSigmaSupport,
                           const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>
                               &igmcorrectionMeiksin = nullptr);
  void resetAsymFitParams();
  void setNominalAmplitude(Float64 ampl) { m_NominalAmplitude = ampl; }

  const std::string &GetName() const { return m_Name; };
  const std::string &GetGroupName() const { return m_GroupName; };
  const Float64 GetNominalAmplitude() const { return m_NominalAmplitude; };
  const std::string &GetVelGroupName() const { return m_VelGroupName; };
  const std::string &GetStrID() const { return m_strID; };

  static EType string2Type(std::string const &s);
  static EForce string2Force(std::string const &s);

protected:
  Int32 m_id = undefIdx;
  EType m_Type = EType::nType_Emission;
  CLineProfile_ptr m_Profile = nullptr;
  EForce m_Force = EForce::nForce_Weak;
  Float64 m_Pos = 0.;
  Float64 m_Offset = 0.;

  std::string m_Name = "";

  // for multiline group
  std::string m_GroupName = "";
  Float64 m_NominalAmplitude = 0.;

  // for offset fitting
  bool m_OffsetFit = false;

  // for velocity fitting
  std::string m_VelGroupName = "";

  std::string m_strID;
};

using CLineVector = std::vector<CLine>;
using CLineMap = std::map<Int32, CLine>;
} // namespace NSEpic

#endif

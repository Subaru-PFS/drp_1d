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
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/ray/ray.h"
#include "RedshiftLibrary/common/exception.h"

using namespace NSEpic;
using namespace std;
#include <fstream>

CRay::CRay(const string& name,
           Float64 pos, UInt32 type,
           CLineProfile_ptr &&profile,
           UInt32 force,
           Float64 amp,
           Float64 width,
           Float64 cut ,
           Float64 posErr,
           Float64 sigmaErr,
           Float64 ampErr,
           const std::string& groupName,
           Float64 nominalAmp,
           const string &velGroupName,
	       Int32 id):
    m_Name(name),
    m_Pos(pos),
    m_Type(type),
    m_Force(force),
    m_Amp(amp),
    m_Width(width),
    m_Cut(cut),
    m_Profile(std::move(profile)),
    m_PosFitErr(posErr),
    m_SigmaFitErr(sigmaErr),
    m_AmpFitErr(ampErr),
    m_GroupName(groupName),
    m_NominalAmplitude(nominalAmp),
    m_VelGroupName(velGroupName),
    m_id(id),
    m_Offset(0.),
    m_OffsetFit(false)
{}

CRay::CRay(const string& name,
           Float64 pos,
	   UInt32 type,
           std::shared_ptr<CLineProfile> profile,
           UInt32 force,
	   Float64 velocityOffset,
	   bool enableVelocityOffset,
           const std::string& groupName,
           Float64 nominalAmp,
           const string &velGroupName,
	   Int32 id):
  m_Name(name),
  m_Pos(pos),
  m_Type(type),
  m_Force(force),
  m_Amp(-1.0),
  m_Width(-1.0),
  m_Cut(-1.0),
  m_Profile(profile->Clone()),
  m_PosFitErr(-1.0),
  m_SigmaFitErr(-1.0),
  m_AmpFitErr(-1.0),
  m_GroupName(groupName),
  m_NominalAmplitude(nominalAmp),
  m_VelGroupName(velGroupName),
  m_id(id),
  m_Offset(velocityOffset),
  m_OffsetFit(enableVelocityOffset)
{
}

CRay::CRay(const string& name,
           Float64 pos, UInt32 type,
           const CLineProfile &profile,
           UInt32 force,
           Float64 amp,
           Float64 width,
           Float64 cut ,
           Float64 posErr,
           Float64 sigmaErr,
           Float64 ampErr,
           const std::string& groupName,
           Float64 nominalAmp,
           const string &velGroupName,
	       Int32 id):
    CRay(   name,
            pos,
            type,
            profile.Clone(),
            force,
            amp,
            width,
            cut ,
            posErr,
            sigmaErr,
            ampErr,
            groupName,
            nominalAmp,
            velGroupName,
	        id)
{}

CRay::CRay(const CRay & other):
    m_id(other.m_id),
    m_Type(other.m_Type),
    m_Profile(other.m_Profile->Clone()),//deep copy for m_Profile
    m_Force(other.m_Force),
    m_Pos(other.m_Pos),
    m_Offset(other.m_Offset),
    m_Amp(other.m_Amp),
    m_Width(other.m_Width),
    m_Cut(other.m_Cut),
    m_PosFitErr(other.m_PosFitErr),
    m_SigmaFitErr(other.m_SigmaFitErr),
    m_AmpFitErr(other.m_AmpFitErr),
    m_Name(other.m_Name),
    m_GroupName(other.m_GroupName),
    m_NominalAmplitude(other.m_NominalAmplitude),
    m_OffsetFit(other.m_OffsetFit),
    m_VelGroupName(other.m_VelGroupName)
    {}



CRay& CRay::operator=(const CRay& other)
{
    m_id = other.m_id;
    m_Type = other.m_Type;
    m_Profile = other.m_Profile->Clone(); //deep copy for m_Profile
    m_Force = other.m_Force;
    m_Pos = other.m_Pos;
    m_Offset = other.m_Offset;
    m_Amp = other.m_Amp;
    m_Width = other.m_Width;
    m_Cut = other.m_Cut;

    m_PosFitErr = other.m_PosFitErr;
    m_SigmaFitErr = other.m_SigmaFitErr;
    m_AmpFitErr = other.m_AmpFitErr;

    m_Name = other.m_Name;

    m_GroupName = other.m_GroupName;
    m_NominalAmplitude = other.m_NominalAmplitude;

    m_OffsetFit = other.m_OffsetFit;

    m_VelGroupName = other.m_VelGroupName;

    return *this;
}


bool CRay::operator < (const CRay& str) const
{
    if(m_Pos == str.m_Pos){
        return (m_Amp < str.m_Amp);
    }else{
        return (m_Pos < str.m_Pos);
    }
}

bool CRay::operator != (const CRay& str) const
{
    if(m_Pos == str.m_Pos){
        return (m_Amp != str.m_Amp);
    }else{
        return (m_Pos != str.m_Pos);
    }
}

void CRay::SetAsymParams(TAsymParams asymParams)
{
    if(!m_Profile)
        throw GlobalException(INTERNAL_ERROR,"CRay::SetAsymParams: lineprofile is not initialized");
    m_Profile->SetAsymParams(asymParams);
}
void CRay::resetAsymFitParams()
{
    if(!m_Profile)
        throw GlobalException(INTERNAL_ERROR,"CRay::resetAsymParams: lineprofile is not initialized");
    m_Profile->resetAsymFitParams();
}
TAsymParams CRay::GetAsymParams() const
{
    if(!m_Profile)
        throw GlobalException(INTERNAL_ERROR,"CRay::GetAsymParams: lineprofile is not initialized");
    return m_Profile->GetAsymParams();
}

bool CRay::GetIsStrong() const
{
    return m_Force == nForce_Strong;
}

bool CRay::GetIsEmission() const
{
    return m_Type == nType_Emission;
}

Int32 CRay::GetType() const
{
    return m_Type;
}

const CLineProfile & CRay::GetProfile() const
{
    if(!m_Profile)
        throw GlobalException(INTERNAL_ERROR,"Current Ray does not have a profile");
    return *m_Profile;
}

void CRay::SetProfile(CLineProfile_ptr&& profile)
{
    m_Profile = std::move(profile);
}

Int32 CRay::GetForce() const
{
    return m_Force;
}

Float64 CRay::GetPosition() const
{
    return m_Pos;
}

Float64 CRay::GetOffset() const
{
    return m_Offset;
}

bool CRay::SetOffset(Float64 val)
{
    if(m_OffsetFit)
    {
        m_Offset = val;
        return true;
    }else{
        return false;
    }
}

bool CRay::GetOffsetFitEnabled() const
{
    return m_OffsetFit;
}

bool CRay::EnableOffsetFit()
{
    m_OffsetFit = true;
    return true;
}

bool CRay::DisableOffsetFit()
{
    m_OffsetFit = false;
    return true;
}

Float64 CRay::GetAmplitude() const
{
    return m_Amp;
}

Float64 CRay::GetWidth() const
{
    return m_Width;
}

Float64 CRay::GetCut() const
{
    return m_Cut;
}

Float64 CRay::GetPosFitError() const
{
    return m_PosFitErr;
}

Float64 CRay::GetSigmaFitError() const
{
    return m_SigmaFitErr;
}

Float64 CRay::GetAmpFitError() const
{
    return m_AmpFitErr;
}

const std::string& CRay::GetName() const
{
    return m_Name;
}

const std::string& CRay::GetGroupName() const
{
    return m_GroupName;
}

const Float64 CRay::GetNominalAmplitude() const
{
    return m_NominalAmplitude;
}

const std::string& CRay::GetVelGroupName() const
{
    return m_VelGroupName;
}

Int32 CRay::GetID() const
{
    return m_id;
}

void CRay::Save(  std::ostream& stream ) const
{
    stream << GetName() << "\t" << GetPosition() << "\t";
    if( GetIsStrong() )
        stream << "Strong"  << "\t";
    else
        stream << "Weak"  << "\t";
    stream << GetCut() << "\t" << GetWidth() << "\t" << GetAmplitude() << "\t" << GetPosFitError() << "\t" << GetSigmaFitError() << "\t" << GetAmpFitError() << "\t";
}

void CRay::SaveDescription(  std::ostream& stream ) const
{
    stream << "#";
    stream << "Name" << "\t" << "Position" << "\t";
    stream << "Force"  << "\t";
    stream << "Cut" << "\t" << "Width" << "\t" << "Amp" << "\t" << "PosFitErr" << "\t" << "SigmaFitErr" << "\t" << "AmpFitErr" << "\t";
}



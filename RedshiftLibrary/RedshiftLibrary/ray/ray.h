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
#ifndef _REDSHIFT_RAY_RAY_
#define _REDSHIFT_RAY_RAY_

#include "RedshiftLibrary/common/datatypes.h"

#include "RedshiftLibrary/ray/lineprofile.h"
#include "RedshiftLibrary/ray/lineprofileLOR.h"
#include "RedshiftLibrary/ray/lineprofileSYM.h"
#include "RedshiftLibrary/ray/lineprofileASYM.h"
#include "RedshiftLibrary/ray/lineprofileASYMFIT.h"

#include <string>
#include <cmath>

namespace NSEpic
{

/**
 * \ingroup Redshift
 * Represent a Single Line.
 */
class CRay
{

public:
    enum EType
    {
        nType_Absorption = 1,
        nType_Emission = 2,
	    nType_All = 3,
    };

    enum EForce
    {
        nForce_Weak = 1,
        nForce_Strong = 2,
    };

    CRay();
    CRay( const std::string& name,
          Float64 pos, UInt32 type,
          std::shared_ptr<CLineProfile> profile,
          UInt32 force,
          Float64 amp=-1.0,
          Float64 width=-1.0,
          Float64 cut=-1.0,
          Float64 posErr=-1.0,
          Float64 sigmaErr=-1.0,
          Float64 ampErr=-1.0,
          const std::string& groupName="-1",
          Float64 nominalAmp=1.0,
          const std::string& velGroupName="-1",
	      Int32 id=-1);

    CRay clone() const;

    bool operator < (const CRay& str) const;
    bool operator != (const CRay& str) const;

    Int32               GetID() const;
    bool                GetIsStrong() const;
    bool                GetIsEmission() const;
    Int32               GetForce() const;
    Int32               GetType() const;
    std::shared_ptr<CLineProfile>        GetProfile() const;
    bool                SetProfile(const std::shared_ptr<CLineProfile>& profile);

    Float64             GetPosition() const;
    Float64             GetOffset() const;
    bool                SetOffset(Float64 val);
    bool                GetOffsetFitEnabled() const;
    bool                EnableOffsetFit();
    bool                DisableOffsetFit();
    
    Float64             GetAmplitude() const;
    Float64             GetWidth() const;
    Float64             GetCut() const;
    Float64             GetPosFitError() const;
    Float64             GetSigmaFitError() const;
    Float64             GetAmpFitError() const;
    const TAsymParams   GetAsymParams();
    void                SetAsymParams(TAsymParams asymParams);
    void                resetAsymFitParams();

    const std::string&  GetName() const;
    const std::string&  GetGroupName() const;
    const Float64 GetNominalAmplitude() const;

    const std::string&  GetVelGroupName() const;

    void                Save( std::ostream& stream ) const;
    void                SaveDescription( std::ostream& stream ) const;

private:
    Int32           m_id = -1;
    Int32           m_Type = 0;
    std::shared_ptr<CLineProfile>    m_Profile=nullptr;
    Int32           m_Force = 0;
    Float64         m_Pos = 0.;
    Float64         m_Offset = 0.;
    Float64         m_Amp = 0.;
    Float64         m_Width = 0.;
    Float64         m_Cut = 0.;

    //fit err
    Float64         m_PosFitErr = 0.;
    Float64         m_SigmaFitErr = 0.;
    Float64         m_AmpFitErr = 0.;

    std::string     m_Name = "";

    //for multiline group
    std::string     m_GroupName = "";
    Float64         m_NominalAmplitude = 0.;

    //for offset fitting
    bool            m_OffsetFit = false;

    //for velocity fitting
    std::string     m_VelGroupName = "";

};

}

#endif

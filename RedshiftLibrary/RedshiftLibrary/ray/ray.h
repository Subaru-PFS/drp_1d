#ifndef _REDSHIFT_RAY_RAY_
#define _REDSHIFT_RAY_RAY_

#include <RedshiftLibrary/common/datatypes.h>

#include <string>
#include <cmath>

namespace NSEpic
{

/**
 * struct that holds ASYMFIXED profile parameters
 */
typedef struct {
    Float64 width, alpha, delta;
} TAsymParams;

/**
 * /ingroup Redshift
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

    enum TProfile
    {
      NONE,
      SYM,
      SYMXL,
      LOR,
      ASYM,
      ASYM2,
      ASYMFIT,
      ASYMFIXED,
      EXTINCT
    };

    typedef std::vector<TProfile> TProfileList;

    CRay();
    CRay( const std::string& name,
          Float64 pos, UInt32 type,
          TProfile profile,
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
          TAsymParams asymParams={NAN, NAN, NAN});

    ~CRay();
    bool operator < (const CRay& str) const;
    bool operator != (const CRay& str) const;

    Bool                GetIsStrong() const;
    Bool                GetIsEmission() const;
    Int32               GetForce() const;
    Int32               GetType() const;
    CRay::TProfile      GetProfile() const;
    bool                SetProfile(CRay::TProfile profile);

    Float64             GetPosition() const;
    Float64             GetOffset() const;
    bool                SetOffset(Float64 val);
    bool                GetOffsetFitEnabled() const;
    bool                EnableOffsetFit(bool val);

    Float64             GetAmplitude() const;
    Float64             GetWidth() const;
    Float64             GetCut() const;
    Float64             GetPosFitError() const;
    Float64             GetSigmaFitError() const;
    Float64             GetAmpFitError() const;
    TAsymParams         GetAsymParams() { return m_asymParams; };
    void                SetAsymParams(TAsymParams asymParams) { m_asymParams = asymParams; };


    const std::string&  GetName() const;
    const std::string&  GetGroupName() const;
    const Float64 GetNominalAmplitude() const;

    const std::string&  GetVelGroupName() const;

    void                Save( std::ostream& stream ) const;
    void                SaveDescription( std::ostream& stream ) const;

    void                ConvertVacuumToAir();

private:

    Int32           m_Type = 0;
    TProfile        m_Profile = NONE;
    Int32           m_Force = 0;
    Float64         m_Pos = 0;
    Float64         m_Offset = 0;
    Float64         m_Amp = 0;
    Float64         m_Width = 0;
    Float64         m_Cut = 0;

    TAsymParams     m_asymParams = {NAN, NAN, NAN};

    //fit err
    Float64         m_PosFitErr = 0;
    Float64         m_SigmaFitErr = 0;
    Float64         m_AmpFitErr = 0;

    std::string     m_Name = "";

    //for multiline group
    std::string     m_GroupName = "";
    Float64         m_NominalAmplitude = 0;

    //for offset fitting
    bool            m_OffsetFit = false;

    //for velocity fitting
    std::string     m_VelGroupName = "";

};


}

#endif

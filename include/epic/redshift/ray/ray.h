#ifndef _REDSHIFT_RAY_RAY_
#define _REDSHIFT_RAY_RAY_

#include <epic/core/common/datatypes.h>

#include <string>

namespace NSEpic
{

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
    };

    enum EForce
    {
        nForce_Weak = 1,
        nForce_Strong = 2,
    };

    CRay();
    CRay( const std::string& name, Float64 pos, UInt32 type, UInt32 force, Float64 amp=-1.0, Float64 width=-1.0, Float64 cut=-1.0, Float64 posErr=-1.0 );
    ~CRay();
    bool operator < (const CRay& str) const;
    bool operator != (const CRay& str) const;

    Bool                GetIsStrong() const;
    Bool                GetIsEmission() const;
    Int32               GetForce() const;
    Int32               GetType() const;
    Float64             GetPosition() const;
    Float64             GetAmplitude() const;
    Float64             GetWidth() const;
    Float64             GetCut() const;
    Float64             GetPosFitError() const;


    const std::string&  GetName() const;
    Void                Save( std::ostream& stream ) const;

    Void                ConvertVacuumToAir();

private:

    Int32           m_Type;
    Int32           m_Force;
    Float64         m_Pos;
    Float64         m_Amp;
    Float64         m_Width;
    Float64         m_Cut;

    //fit err
    Float64         m_PosFitErr;

    std::string     m_Name;
};


}

#endif

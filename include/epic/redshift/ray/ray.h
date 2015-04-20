#ifndef _REDSHIFT_RAY_RAY_
#define _REDSHIFT_RAY_RAY_

#include <epic/core/common/datatypes.h>

#include <string>

namespace NSEpic
{

/**
 * Represent a Single Ray.
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
    CRay( const std::string& name, Float64 pos, UInt32 type );
    ~CRay();

    Bool                GetIsStrong() const;
    Bool                GetIsEmission() const;
    Float64             GetPosition() const;
    const std::string&  GetName() const;
    const std::string GetDescription() const;

private:

    UInt32         m_Type;
    UInt32         m_Force;
    Float64        m_Pos;
    std::string    m_Name;
};


}

#endif

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
        nType_Absorption = (1<<0),
        nType_Emission = (1<<1),
        nType_Weak = (1<<2),
        nType_Strong = (1<<3),
    };

    CRay();
    CRay( const std::string& name, Float64 pos, UInt32 type );
    ~CRay();

    Bool                GetIsStrong() const;
    Bool                GetIsEmission() const;
    Float64             GetPosition() const;
    const std::string&  GetName() const;

private:

    UInt32         m_Type;
    Float64        m_Pos;
    std::string    m_Name;
};


}

#endif

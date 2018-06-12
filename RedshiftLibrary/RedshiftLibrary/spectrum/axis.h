#ifndef _REDSHIFT_SPECTRUM_AXIS_
#define _REDSHIFT_SPECTRUM_AXIS_

#include <RedshiftLibrary/common/datatypes.h>


#include <vector>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CSpectrumAxis
{

public:

    CSpectrumAxis();
    CSpectrumAxis( UInt32 n );
    CSpectrumAxis( const Float64* samples, UInt32 n );
    ~CSpectrumAxis();

    CSpectrumAxis& operator=(const CSpectrumAxis& other);

    Float64 operator[]( const UInt32 i ) const;
    Float64& operator[]( const UInt32 i );

    const Float64*      GetSamples() const;
    Float64*            GetSamples();
    UInt32              GetSamplesCount() const;

    virtual Void        SetSize( UInt32 s );

protected:

    TAxisSampleList     m_Samples;

};

inline
Float64 CSpectrumAxis::operator[]( const UInt32 i ) const
{
    return m_Samples[i];
}

inline
Float64& CSpectrumAxis::operator[]( const UInt32 i )
{
    return m_Samples[i];
}

inline
UInt32 CSpectrumAxis::GetSamplesCount() const
{
    return m_Samples.size();
}


inline
Float64* CSpectrumAxis::GetSamples()
{
    return m_Samples.data();
}

inline
const Float64* CSpectrumAxis::GetSamples() const
{
    return m_Samples.data();
}


}

#endif

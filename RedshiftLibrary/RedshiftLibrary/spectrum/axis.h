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
    explicit CSpectrumAxis( UInt32 n );
    CSpectrumAxis( const Float64* samples, UInt32 n );
    ~CSpectrumAxis();

    CSpectrumAxis& operator=(const CSpectrumAxis& other);
    CSpectrumAxis& operator*=(Float64 op);
    Float64 operator[]( const UInt32 i ) const;
    Float64& operator[]( const UInt32 i );

    const Float64*           GetSamples() const;
    Float64*                 GetSamples();
    const TAxisSampleList&   GetSamplesVector() const;
    TAxisSampleList&         GetSamplesVector();
    UInt32                   GetSamplesCount() const;
    UInt32                   GetSamplesCount();
    virtual void        SetSize( UInt32 s );

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

inline
TAxisSampleList& CSpectrumAxis::GetSamplesVector()
{
    return m_Samples;
}

inline
const TAxisSampleList& CSpectrumAxis::GetSamplesVector() const
{
    return m_Samples;
}

}
#endif

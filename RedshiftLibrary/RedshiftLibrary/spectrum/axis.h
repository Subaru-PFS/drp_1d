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

    CSpectrumAxis() = default;
    CSpectrumAxis(const CSpectrumAxis & other) = default;
    CSpectrumAxis(CSpectrumAxis && other) = default;
    explicit CSpectrumAxis( UInt32 n, Float64 value = 0.0 );
    CSpectrumAxis( const Float64* samples, UInt32 n );
    ~CSpectrumAxis() = default;
    CSpectrumAxis& operator=(const CSpectrumAxis& other) = default;
    CSpectrumAxis& operator=(CSpectrumAxis&& other) = default;
    CSpectrumAxis& operator*=(const Float64 op);
    Float64 operator[]( const UInt32 i ) const;
    Float64& operator[]( const UInt32 i );
    
    Int32 MaskAxis(TFloat64List& mask, CSpectrumAxis& maskedAxis) const;

    const Float64*           GetSamples() const;
    Float64*                 GetSamples();
    const TAxisSampleList&   GetSamplesVector() const;
    TAxisSampleList&         GetSamplesVector();
    UInt32                   GetSamplesCount() const;
    UInt32                   GetSamplesCount();
    virtual void             SetSize( UInt32 s );
    Int32 extractFrom(const CSpectrumAxis& other, Int32 startIdx, Int32 endIdx);
    Bool isEmpty() const ;
protected:

    TAxisSampleList          m_Samples;

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
inline
Bool CSpectrumAxis::isEmpty() const{
    return m_Samples.size()==0;
}
}
#endif

#ifndef _REDSHIFT_SPECTRUM_LSF_
#define _REDSHIFT_SPECTRUM_LSF_

#include <RedshiftLibrary/common/datatypes.h>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CLSF
{

public:

    CLSF(const Float64 sigma=0.0);
    ~CLSF();

    const Float64       GetSigma() const;
    void                SetSigma(const Float64 sigma);
    bool                IsValid() const;

private:

    Float64             m_sigma=0.0;
};

}

#endif

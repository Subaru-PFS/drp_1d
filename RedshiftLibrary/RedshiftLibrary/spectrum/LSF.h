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

    virtual ~CLSF(){};

    virtual const Float64       GetSigma() const=0;
    virtual void                SetSigma(const Float64 sigma)=0;
    virtual bool                IsValid() const=0;

};

}

#endif

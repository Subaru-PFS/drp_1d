#ifndef _REDSHIFT_COMMON_GAUSSIANFIT_
#define _REDSHIFT_COMMON_GAUSSIANFIT_

#include <epic/redshift/common/datatypes.h>

namespace __NS__
{

class CGaussianFit
{

public:

    CGaussianFit( );
    ~CGaussianFit();

    Bool    Compute( const Float64* x, const Float64* y, Int32 n, Int32 polyOrder );

private:

};

}

#endif

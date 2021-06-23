#ifndef _REDSHIFT_OPERATOR_SPECTRARESULT_
#define _REDSHIFT_OPERATOR_SPECTRARESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>


using namespace std;
namespace NSEpic
{

class CSpectraFluxResult : public COperatorResult
{

public:

    CSpectraFluxResult();
    virtual ~CSpectraFluxResult();

    TFloat64List   fluxes;
    TFloat64List   wavel;


};


}

#endif

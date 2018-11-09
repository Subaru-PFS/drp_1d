#ifndef _REDSHIFT_OPERATOR_CONTINUUMMODELSOLUTION_
#define _REDSHIFT_OPERATOR_CONTINUUMMODELSOLUTION_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/continuum/indexes.h>

namespace NSEpic
{


struct CContinuumModelSolution
{

    //template continuum
    std::string         tplName;
    Float64             tplAmplitude;
    Float64             tplMerit;
    Float64             tplDustCoeff;
    Int32               tplMeiksinIdx;

    //polynom
    Int32           pOrder;
    TFloat64List    pCoeffs;
};


}

#endif

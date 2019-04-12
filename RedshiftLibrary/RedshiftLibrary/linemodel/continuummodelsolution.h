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
    Float64             tplDustCoeff;
    Int32               tplMeiksinIdx;
    Float64             tplRedshift;

    Float64             tplMerit;
    Float64             tplDtm;
    Float64             tplMtm;
    Float64             tplLogPrior;

    //polynom
    TFloat64List    pCoeffs;
};


}

#endif

#ifndef _REDSHIFT_LINEMODEL_CONTINUUMMODELSOLUTION_
#define _REDSHIFT_LINEMODEL_CONTINUUMMODELSOLUTION_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/continuum/indexes.h>
#include <cmath>

namespace NSEpic
{


struct CContinuumModelSolution
{

    //template continuum
    std::string         tplName;
    Float64             tplAmplitude = NAN;
    Float64             tplAmplitudeError = NAN;
    Float64             tplDustCoeff = NAN;
    Int32               tplMeiksinIdx = -1;
    Float64             tplRedshift = NAN;

    Float64             tplMerit = NAN;
    Float64             tplDtm = NAN;
    Float64             tplMtm = NAN;
    Float64             tplLogPrior = NAN;

    //polynom
    TFloat64List    pCoeffs;
};


}

#endif

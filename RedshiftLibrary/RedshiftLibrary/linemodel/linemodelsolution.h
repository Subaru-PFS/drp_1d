#ifndef _REDSHIFT_LINEMODEL_LINEMODELSOLUTION_
#define _REDSHIFT_LINEMODEL_LINEMODELSOLUTION_

#include "RedshiftLibrary/processflow/result.h"
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/continuum/indexes.h"

namespace NSEpic
{


struct CLineModelSolution
{

    std::vector<Float64> ElementId;     //id of the linemodel element it is part of
    std::vector<Float64> Amplitudes;
    std::vector<CRay> Rays;
    std::vector<Float64> Errors;    //noise sigma
    std::vector<Float64> FittingError;    //ModelLeastSquare error under each line
    std::vector<Float64> CenterContinuumFlux;    //Continuum flux value at the center of each line
    std::vector<Float64> ContinuumError;    //Continuum error value for each line
    std::vector<Float64> Sigmas;    //width for each line
    std::vector<Float64> Fluxs;    //Flux for each line
    std::vector<Float64> FluxErrors;    //Flux error for each line
    std::vector<Float64> FluxDirectIntegration;    //Flux obtained by direct integration for each line

    Float64 snrHa;
    Float64 lfHa;
    Float64 snrOII;
    Float64 lfOII;
    Int32 NLinesAboveSnrCut;

    std::vector<Float64> LambdaObs;  //observed position in Angstrom
    std::vector<Float64> Velocity;  //dispersion velocity in km/s
    std::vector<Float64> Offset;    //line offset in km/s
    std::vector<Bool> OutsideLambdaRange;
    std::vector<TInt32Range> fittingIndexRange;
    std::vector<std::string> fittingGroupInfo;

    Float64 LyaWidthCoeff;
    Float64 LyaAlpha;
    Float64 LyaDelta;

    Float64 AbsorptionVelocity;
    Float64 EmissionVelocity;
    Float64 Redshift;

    Int32 nDDL;

};


}

#endif

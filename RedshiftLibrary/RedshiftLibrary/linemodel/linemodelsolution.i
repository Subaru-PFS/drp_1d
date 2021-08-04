
  class CLineModelSolution: public COperatorResult
{
 public:

  CLineModelSolution();
 
  
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
    std::vector<Float64> FluxDirectIntegrationError;    //Flux obtained by direct integration for each line
    std::vector<Int32> rayId;

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

    void fillRayIds();
      
};

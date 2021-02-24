class CModelFittingResult : public COperatorResult
{

public:

    CModelFittingResult(CLineModelSolution _lineModelSolution, Float64 _redshift, Float64 _merit, CRayCatalog::TRayVector _restRayList, Float64 _velEmission=-1.0, Float64 _velAbsorption=-1.0 );
    CModelFittingResult();
    virtual ~CModelFittingResult();

  void Load( const char* filePath );

    const CLineModelSolution& GetLineModelSolution() const;

  //copies from restRayList For output only
    std::vector<Int32> rayId;
    std::vector<Float64> rayLambdaRest;
    std::vector<Float64> FittedRaysFlux;
    std::vector<Float64> FittedRaysLambda;

private:

    CLineModelSolution LineModelSolution;
    Float64 Redshift;
    Float64 Merit;

    CRayCatalog::TRayVector restRayList;
    Float64 VelocityEmission;
    Float64 VelocityAbsorption;


};

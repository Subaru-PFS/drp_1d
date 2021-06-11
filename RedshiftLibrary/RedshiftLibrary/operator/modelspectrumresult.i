  class CModelSpectrumResult : public COperatorResult
{

 public:

  CModelSpectrumResult(const CSpectrum& spc);
  CModelSpectrumResult(CSpectrum&& spc);

  //  CModelSpectrumResult();
    
   const TFloat64List ModelLambda;
   const TFloat64List ModelFlux;

};

  class CModelSpectrumResult : public COperatorResult
{

public:

    CModelSpectrumResult(const CSpectrum& spc);
    //  CModelSpectrumResult();
    virtual ~CModelSpectrumResult();

    CSpectrum& GetSpectrum();

    const TFloat64List& ModelLambda;
    const TFloat64List& ModelFlux;
private:
  CSpectrum m_model;
};

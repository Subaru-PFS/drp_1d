#ifndef _REDSHIFT_HYBRID_ELEMENTLIST_
#define _REDSHIFT_HYBRID_ELEMENTLIST_

#include "RedshiftLibrary/linemodel/abstractfitter.h"

namespace NSEpic

{

// class CRegulament;
class CHybridFitter : public CAbstractFitter {
public:
  CHybridFitter(CLineModelElementList &elements,
                std::shared_ptr<const CSpectrum> inputSpectrum,
                std::shared_ptr<const TLambdaRange> lambdaRange,
                std::shared_ptr<CSpectrumModel> spectrumModel);

  void fit(Float64 redshift);

  Int32 m_cont_reestim_iterations = 0;

  bool m_opt_enable_improveBalmerFit = false;

private:
  //  CRegulament &m_Regulament;
  Int32 fitAmplitudesHybrid(Float64 redshift);
  Int32 improveBalmerFit(Float64 redshift);
};
} // namespace NSEpic
#endif

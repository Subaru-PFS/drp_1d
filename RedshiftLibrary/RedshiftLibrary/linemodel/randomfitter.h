#ifndef _REDSHIFT_RANDOM_ELEMENTLIST_
#define _REDSHIFT_RANDOM_ELEMENTLIST_

#include "RedshiftLibrary/linemodel/abstractfitter.h"

namespace NSEpic

{

// class CRegulament;
class CRandomFitter : public CAbstractFitter {
public:
  CRandomFitter(CLineModelElementList &elements,
                std::shared_ptr<const CSpectrum> inputSpectrum,
                std::shared_ptr<const TLambdaRange> lambdaRange,
                std::shared_ptr<CSpectrumModel> spectrumModel);

  void fit(Float64 redshift);

private:
  Float64 getContinuumMeanUnderElement(Int32 eltId) const;

  //  CRegulament &m_Regulament;
};
} // namespace NSEpic
#endif

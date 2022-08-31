#ifndef _REDSHIFT_INDIVIDUAL_ELEMENTLIST_
#define _REDSHIFT_INDIVIDUAL_ELEMENTLIST_

#include "RedshiftLibrary/linemodel/abstractfitter.h"

namespace NSEpic

{

// class CRegulament;
class CIndividualFitter : public CAbstractFitter {
public:
  CIndividualFitter(CLineModelElementList &elements,
                    std::shared_ptr<const CSpectrum> inputSpectrum,
                    std::shared_ptr<const TLambdaRange> lambdaRange,
                    std::shared_ptr<CSpectrumModel> spectrumModel);

  void fit(Float64 redshift);

private:
  //  CRegulament &m_Regulament;
};
} // namespace NSEpic
#endif

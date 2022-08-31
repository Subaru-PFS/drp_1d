#ifndef _REDSHIFT_SVD_ELEMENTLIST_
#define _REDSHIFT_SVD_ELEMENTLIST_

#include "RedshiftLibrary/linemodel/abstractfitter.h"

namespace NSEpic

{

class CSvdFitter : public CAbstractFitter {
public:
  CSvdFitter(CLineModelElementList &elements,
             std::shared_ptr<const CSpectrum> inputSpectrum,
             std::shared_ptr<const TLambdaRange> lambdaRange,
             std::shared_ptr<CSpectrumModel> spectrumModel);

  void fit(Float64 redshift);
};
} // namespace NSEpic
#endif

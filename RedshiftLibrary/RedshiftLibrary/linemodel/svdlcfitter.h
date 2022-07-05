#ifndef _REDSHIFT_SVDLC_ELEMENTLIST_
#define _REDSHIFT_SVDLC_ELEMENTLIST_

#include "RedshiftLibrary/linemodel/abstractfitter.h"

namespace NSEpic

{
class CTemplateCatalog;
// class CRegulament;
class CSvdlcFitter : public CAbstractFitter {
public:
  CSvdlcFitter(CLineModelElementList &elements,
               std::shared_ptr<const CSpectrum> inputSpectrum,
               std::shared_ptr<const TLambdaRange> lambdaRange,
               std::shared_ptr<CSpectrumModel> spectrumModel,
               Int32 polyOrder = -1);

  void fit(Float64 redshift);

  Int32 m_fitc_polyOrder = -1;
  std::shared_ptr<const CTemplateCatalog> m_tplCatalog;
  TStringList m_tplCategoryList;

private:
  Int32 fitAmplitudesLinesAndContinuumLinSolve(
      const TInt32List &EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
      const CSpectrumFluxAxis &fluxAxis,
      const CSpectrumFluxAxis &continuumfluxAxis, TFloat64List &ampsfitted,
      TFloat64List &errorsfitted, Float64 &chisquare, Float64 redshift,
      Int32 polyOrder = -1);
};
} // namespace NSEpic
#endif

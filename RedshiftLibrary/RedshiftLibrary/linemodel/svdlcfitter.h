// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_SVDLC_ELEMENTLIST_
#define _REDSHIFT_SVDLC_ELEMENTLIST_

#include "RedshiftLibrary/linemodel/abstractfitter.h"
#include "RedshiftLibrary/linemodel/continuummanager.h"

namespace NSEpic

{

// class CRegulament;
class CSvdlcFitter : public CAbstractFitter {
public:
  CSvdlcFitter(const std::shared_ptr<CElementsLists> &elementsVector,
               const CCSpectrumVectorPtr &inputSpcs,
               const CTLambdaRangePtrVector &lambdaRanges,
               const CSpcModelVectorPtr &spectrumModels,
               const CLineMap &restLineList,
               const std::shared_ptr<Int32> &curObsPtr,
               std::shared_ptr<CContinuumManager> continuumManager,
               Int32 polyOrder = -1, bool enableAmplitudeOffsets = false,
               bool enableLambdaOffsetsFit = false);

private:
  void doFit(Float64 redshift) override;
  Int32 fitAmplitudesLinesAndContinuumLinSolve(
      const TInt32List &EltsIdx, const CSpectrumSpectralAxis &spectralAxis,
      TFloat64List &ampsfitted, TFloat64List &errorsfitted, Float64 &chisquare,
      Float64 redshift);

  void fillMatrix(Int32 imin, Int32 imax, Float64 redshift,
                  const TInt32List &EltsIdx,
                  const CSpectrumSpectralAxis &spectralAxis,
                  const CSpectrumFluxAxis &continuumfluxAxis,
                  gsl_matrix *X) const;
  Int32 m_fitc_polyOrder = -1;

  std::shared_ptr<CContinuumManager> m_continuumManager;
  const CSpectrumSpectralAxis &m_spectralAxis;
};
} // namespace NSEpic
#endif

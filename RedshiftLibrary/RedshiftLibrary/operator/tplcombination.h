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
#ifndef _REDSHIFT_OPERATOR_TPLCOMBINATION_
#define _REDSHIFT_OPERATOR_TPLCOMBINATION_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/operator/operator.h"
#include "RedshiftLibrary/operator/templatefitting.h"
#include "RedshiftLibrary/processflow/resultstore.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/statistics/priorhelper.h"

#include <gsl/gsl_matrix_double.h>

namespace NSEpic {
class CSpectrum;
class COperatorResult;
class CModelSpectrumResult;

struct STplcombination_basicfitresult : TFittingIsmIgmResult {
  STplcombination_basicfitresult(Int32 EbmvListSize, Int32 MeiksinListSize,
                                 Int32 componentCount)
      : TFittingIsmIgmResult(EbmvListSize, MeiksinListSize),
        fittingAmplitudes(componentCount, NAN),
        fittingAmplitudeErrors(componentCount, NAN),
        fittingAmplitudeSigmas(componentCount, NAN),
        fittingAmplitudesInterm(
            EbmvListSize,
            std::vector<TFloat64List>(MeiksinListSize,
                                      TFloat64List(componentCount, NAN))),
        tplNames(componentCount),
        COV(componentCount, TFloat64List(componentCount, NAN)){};

  TFloat64List fittingAmplitudes;
  TFloat64List fittingAmplitudeErrors;
  TFloat64List fittingAmplitudeSigmas;

  std::vector<std::vector<TFloat64List>>
      fittingAmplitudesInterm; // intermediate amplitudes
  TStringList tplNames;        // cause combination of templates

  Float64 SNR = NAN;
  std::vector<TFloat64List> COV;
};

class COperatorTplcombination {
public:
  std::shared_ptr<COperatorResult>
  Compute(const CSpectrum &spectrum, const TTemplateConstRefList &tplList,
          const TFloat64Range &lambdaRange, const TFloat64List &redshifts,
          Float64 overlapThreshold,
          const std::vector<CMask> &additional_spcMasks,
          const std::string &opt_interp, bool opt_extinction = false,
          bool opt_dustFitting = false,
          const CPriorHelper::TPriorZEList &logpriorze =
              CPriorHelper::TPriorZEList(),
          Int32 FitEbmvIdx = undefIdx, Int32 FitMeiksinIdx = undefIdx);

  Float64 ComputeDtD(const CSpectrumFluxAxis &spcFluxAxis,
                     const TInt32Range &range); // could be also made static
  std::shared_ptr<CModelSpectrumResult> ComputeSpectrumModel(
      const CSpectrum &spectrum, const TTemplateConstRefList &tplList,
      Float64 redshift, Float64 EbmvCoeff, Int32 meiksinIdx,
      const TFloat64List &amplitudes, const TFloat64Range &lambdaRange,
      const Float64 overlapThreshold);

private:
  void BasicFit_preallocateBuffers(const CSpectrum &spectrum,
                                   const TTemplateConstRefList &tplList);

  void BasicFit(const CSpectrum &spectrum, const TTemplateConstRefList &tplList,
                const TFloat64Range &lambdaRange, Float64 redshift,
                Float64 overlapThreshold,
                STplcombination_basicfitresult &fittingResults,
                Float64 forcedAmplitude, bool opt_extinction,
                bool opt_dustFitting, CMask spcMaskAdditional,
                const CPriorHelper::TPriorEList &logpriore,
                const TInt32List &MeiksinList, const TInt32List &EbmvList);
  void RebinTemplate(const CSpectrum &spectrum,
                     const TTemplateConstRefList &tplList, Float64 redshift,
                     const TFloat64Range &lambdaRange,
                     TFloat64Range &currentRange, Float64 &overlapRate,
                     const Float64 overlapThreshold);
  // buffers for the interpolated axis (templates & spectrum)
  std::vector<CTemplate> m_templatesRebined_bf;
  std::vector<CMask> m_masksRebined_bf;
  CSpectrumSpectralAxis m_spcSpectralAxis_restframe;

  Float64 EstimateLikelihoodCstLog(const CSpectrum &spectrum,
                                   const TFloat64Range &lambdaRange);

  Float64 ComputeXi2_bruteForce(const CSpectrumFluxAxis &correctedFlux,
                                const CSpectrumFluxAxis &spcFluxAxis,
                                const Int32 imin_lbda);
  Float64 GetNormFactor(const CSpectrumFluxAxis spcFluxAxis, Int32 kStart,
                        Int32 n);
};

} // namespace NSEpic

#endif

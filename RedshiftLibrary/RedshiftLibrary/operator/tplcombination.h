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
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/operator.h"

#include "RedshiftLibrary/operator/modelspectrumresult.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/processflow/resultstore.h"

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/statistics/priorhelper.h"
#include <gsl/gsl_matrix_double.h>
namespace NSEpic
{
class CSpectrum;
class COperatorResult;
class CModelSpectrumResult;

class COperatorTplcombination
{
public:

    std::shared_ptr<COperatorResult> Compute(const CSpectrum& spectrum,
                                             const TTemplateConstRefList& tplList,
                                             const TFloat64Range& lambdaRange,
                                             const TFloat64List& redshifts,
                                             Float64 overlapThreshold,
                                             std::vector<CMask> additional_spcMasks, 
                                             std::string opt_interp, 
                                             Int32 opt_extinction=0, 
                                             Int32 opt_dustFitting=0,
                                             CPriorHelper::TPriorZEList logpriorze=CPriorHelper::TPriorZEList(),
                                             Bool keepigmism = false,
                                             Float64 FitEbmvCoeff=-1.,
                                             Int32 FitMeiksinIdx=-1);

    Float64 ComputeDtD(const CSpectrumFluxAxis& spcFluxAxis, const TInt32Range& range); //could be also made static
    std::shared_ptr<CModelSpectrumResult>   ComputeSpectrumModel(   
                                    const CSpectrum& spectrum,
                                    const TTemplateConstRefList& tplList,
                                    Float64 redshift,
                                    Float64 EbmvCoeff,
                                    Int32 meiksinIdx,
                                    const TFloat64List& amplitudes,
                                    std::string opt_interp,
                                    const TFloat64Range& lambdaRange,
                                    const Float64 overlapThreshold);
private:

    struct STplcombination_basicfitresult
    {
        COperator::EStatus status;
        Float64     overlapRate;
        Float64     chisquare;
        TFloat64List    fittingAmplitudes;
        TFloat64List    fittingAmplitudeErrors;
        TFloat64List    fittingAmplitudeSigmas;
        std::vector<TFloat64List>    ChiSquareInterm;
        std::vector<TFloat64List>    IsmCalzettiCoeffInterm;
        std::vector<TInt32List>      IgmMeiksinIdxInterm;
        std::vector<std::vector<TFloat64List>>    fittingAmplitudesInterm; //intermediate amplitudes
        std::vector<std::string> tplNames; //cause combination of templates
        Int32 IGMIdx;
        Float64 EbmvCoeff;
        Float64 SNR;
        std::vector<TFloat64List> COV; 
    };

    void BasicFit_preallocateBuffers(const CSpectrum& spectrum,
                                     const TTemplateConstRefList& tplList);

    void BasicFit(const CSpectrum& spectrum,
                  const TTemplateConstRefList& tplList,
                  const TFloat64Range& lambdaRange,
                  Float64 redshift,
                  Float64 overlapThreshold,
                  STplcombination_basicfitresult& fittingResults,
                  std::string opt_interp, Float64 forcedAmplitude=-1, 
                  Int32 opt_extinction=0, 
                  Int32 opt_dustFitting=0, 
                  CMask spcMaskAdditional=CMask(),
                  CPriorHelper::TPriorEList logpriore=CPriorHelper::TPriorEList(),
                  const TInt32List& MeiksinList=TInt32List(-1),
                  const TInt32List& EbmvList=TInt32List(-1));
    void RebinTemplate(const CSpectrum& spectrum,
                        const TTemplateConstRefList& tplList,
                        Float64 redshift,
                        const TFloat64Range& lambdaRange,
                        std::string opt_interp,
                        TFloat64Range& currentRange,
                        Float64& overlapRate,
                        const Float64 overlapThreshold);
    // buffers for the interpolated axis (templates & spectrum)
    std::vector<CTemplate>   m_templatesRebined_bf;
    std::vector<CMask>       m_masksRebined_bf;
    CSpectrumSpectralAxis    m_spcSpectralAxis_restframe;

    Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);

    Float64 ComputeXi2_bruteForce(const CSpectrumFluxAxis& correctedFlux, 
                                  const CSpectrumFluxAxis& spcFluxAxis,
                                  const Int32 imin_lbda);
    Float64 GetNormFactor(const CSpectrumFluxAxis spcFluxAxis, UInt32 kStart, UInt32 n);
};


}

#endif

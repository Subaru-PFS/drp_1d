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
#ifndef _REDSHIFT_OPERATOR_TEMPLATE_FITTING_
#define _REDSHIFT_OPERATOR_TEMPLATE_FITTING_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/operator/templatefittingBase.h"
#include "RedshiftLibrary/operator/templatefittingresult.h"
#include "RedshiftLibrary/common/mask.h"

#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/template/template.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/statistics/priorhelper.h"
#include "RedshiftLibrary/operator/modelspectrumresult.h"

#include <numeric>
namespace NSEpic
{
struct TFittingResult {
    Float64 chiSquare = INFINITY;
    Float64 chiSquare_phot = 0.0; // should be initialized at zero to be summed unconditionnaly 
    Float64 ampl = NAN;
    Float64 ampl_err = NAN;
    Float64 ampl_sigma = NAN;
    Float64 sumCross = 0.;
    Float64 sumT = 0.;
    Float64 sumS = 0.;
    Float64 logprior = 0.;
};

struct TFittingIsmIgmResult : TFittingResult{
    TFittingIsmIgmResult(Int32 EbmvListSize, Int32 MeiksinListSize):
            ChiSquareInterm(EbmvListSize, TFloat64List(MeiksinListSize, DBL_MAX)),
            IsmCalzettiCoeffInterm(EbmvListSize, TFloat64List(MeiksinListSize, NAN)),
            IgmMeiksinIdxInterm(EbmvListSize, TInt32List(MeiksinListSize, -1)) {}

    Float64 overlapRate = NAN;
    Float64 EbmvCoeff = NAN;
    Int32 MeiksinIdx = -1;
    std::vector<TFloat64List> ChiSquareInterm;
    std::vector<TFloat64List> IsmCalzettiCoeffInterm;
    std::vector<TInt32List> IgmMeiksinIdxInterm;
    COperator::EStatus status = COperator::EStatus::nStatus_DataError;
};

class COperatorTemplateFitting : public COperatorTemplateFittingBase
{

public:

    COperatorTemplateFitting(const CSpectrum & spectrum, const TFloat64Range& lambdaRange, const TFloat64List & redshifts=TFloat64List()):
        COperatorTemplateFittingBase(spectrum, lambdaRange, redshifts) {};
    COperatorTemplateFitting(const CSpectrum && spectrum, const TFloat64Range& lambdaRange, const TFloat64List & redshifts) = delete;
    virtual ~COperatorTemplateFitting() = default;
    
    std::shared_ptr<COperatorResult> Compute( const std::shared_ptr<const CTemplate> & tpl,
                                              Float64 overlapThreshold,
                                              const std::vector<CMask> &additional_spcMasks,
                                              std::string opt_interp,
                                              Int32 opt_extinction=0,
                                              Int32 opt_dustFitting=-1,
                                              const CPriorHelper::TPriorZEList &logpriorze=CPriorHelper::TPriorZEList(),
                                              bool keepigmism = false,
                                              Float64 FitEbmvCoeff=-1.,
                                              Int32 FitMeiksinIdx=-1) override;

protected:

    TFittingIsmIgmResult BasicFit(const std::shared_ptr<const CTemplate>& tpl,
                  Float64 redshift,
                  Float64 overlapThreshold,
                  std::string opt_interp,
                  Float64 forcedAmplitude,
                  Int32 opt_extinction,
                  Int32 opt_dustFitting,
                  CMask spcMaskAdditional,
                  const CPriorHelper::TPriorEList & logpriore,
                  const TInt32List& MeiksinList,
                  const TInt32List& EbmvList);

    virtual void InitIsmIgmConfig( Float64 redshift,
                           const std::shared_ptr<const CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti,
                           const std::shared_ptr<const CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin,
                           Int32 EbmvListSize);

    virtual bool CheckLyaIsInCurrentRange(const TFloat64Range & currentRange) const {
        return  currentRange.GetBegin() > 1216.0;
    };

    virtual bool ApplyMeiksinCoeff(Int32 meiksinIdx) {
        return m_templateRebined_bf.ApplyMeiksinCoeff(meiksinIdx);};

    virtual bool ApplyDustCoeff(Int32 kEbmv) {
        return m_templateRebined_bf.ApplyDustCoeff(kEbmv);};

    virtual TFittingResult ComputeLeastSquare(  Int32 kM,
                                                Int32 kEbmv,
                                                const CPriorHelper::SPriorTZE & logprior,
                                                const CMask & spcMaskAdditional);
    
    TFittingResult ComputeCrossProducts(Int32 kM,
                                        Int32 kEbmv_,
                                        const CMask & spcMaskAdditional);

    void ComputeAmplitudeAndChi2( TFittingResult & fitres, const CPriorHelper::SPriorTZE & logpriorTZ) const;

    bool m_option_igmFastProcessing;
    Int32 m_kStart, m_kEnd;
    Float64 m_forcedAmplitude;

    TFloat64List m_sumCross_outsideIGM;
    TFloat64List m_sumT_outsideIGM;
    TFloat64List m_sumS_outsideIGM;

    bool m_amplForcePositive=true;

};


}

#endif

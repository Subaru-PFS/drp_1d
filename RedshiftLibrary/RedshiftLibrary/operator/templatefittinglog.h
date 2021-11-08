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
#ifndef _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_
#define _REDSHIFT_OPERATOR_CHISQUARELOGLAMBDA_

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

#include <fftw3.h>

namespace NSEpic
{

class COperatorTemplateFittingLog : public COperatorTemplateFittingBase
{

public:
    COperatorTemplateFittingLog();
    ~COperatorTemplateFittingLog();

    COperatorTemplateFittingLog(const COperatorTemplateFittingLog & other) = delete; 
    COperatorTemplateFittingLog(COperatorTemplateFittingLog && other) = delete; 
    COperatorTemplateFittingLog& operator=(const COperatorTemplateFittingLog& other) = delete;  
    COperatorTemplateFittingLog& operator=(COperatorTemplateFittingLog&& other) = delete; 

    std::shared_ptr<COperatorResult> Compute( const CSpectrum& rebinnedSpectrum,
                                              const CTemplate& rebinnedTpl,
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

    inline  bool IsFFTProcessing() override{return true;}; 

    //made public for unit-testing
    TInt32Range FindTplSpectralIndex(const TFloat64Range & redshiftrange) const;
    TInt32Range FindTplSpectralIndex( const CSpectrumSpectralAxis & spcSpectrailAxis, 
                                      const CSpectrumSpectralAxis& tplSpectralAxis,
                                      const TFloat64Range & redshiftrange) const;
    //log grid data
    CTemplate       m_templateRebinedLog;
    CSpectrum       m_spectrumRebinedLog;

private:
    //hardcoded config: FIT_RANGEZ
    bool verboseLogFitFitRangez = false;
    bool verboseExportFitRangez = false;
    bool verboseExportFitRangez_model = false;
    UInt32 exportIGMIdx = 5;
    UInt32 exportISMIdx = -1;

    //hardcoded config: XTY_FFT
    bool verboseLogXtYFFT = false;
    bool verboseExportXtYFFT = false;

    Int32 FitAllz(const TFloat64Range& lambdaRange,
                  std::shared_ptr<CTemplateFittingResult> result,
                  TInt32List MeiksinList=TInt32List(1, 0),
                  TInt32List EbmvList=TInt32List(1, 0),
                  CMask spcMaskAdditional=CMask(),
                  CPriorHelper::TPriorZEList logpriorze=CPriorHelper::TPriorZEList());

    Int32 FitRangez(const TFloat64List & inv_err2,
                    TInt32Range& range,
                    std::shared_ptr<CTemplateFittingResult> result,
                    TInt32List MeiksinList,
                    TInt32List EbmvList,
                    const Float64& dtd);

    Int32 EstimateXtY(const TFloat64List& X, const TFloat64List& Y,
                      UInt32 nshifts, TFloat64List& XtY, Int32 precomputedFFT=-1);
    Int32 InitFFT(Int32 n);
    Int32 EstimateXtYSlow(const TFloat64List& X, const TFloat64List& Y, UInt32 nShifts,
                          TFloat64List& XtY);
    Int32 EstimateMtMFast(const TFloat64List& X, const TFloat64List& Y, UInt32 nShifts, TFloat64List& XtY);

    Int32 InterpolateResult(const TFloat64List& in, TFloat64List& inGrid,
                            const TFloat64List& tgtGrid, TFloat64List& out, Float64 defaultValue);

    void freeFFTPlans();
    void freeFFTPrecomputedBuffers();

    Bool m_enableISM = 1;
    Bool m_enableIGM = 1; 

    //buffers for fft computation
    Int32 m_nPaddedSamples;
    Float64 *inSpc;
    fftw_complex *outSpc;
    fftw_plan pSpc;
    Float64 *inTpl;
    Float64 *inTpl_padded;
    fftw_complex *outTpl;
    fftw_plan pTpl;
    fftw_complex* outCombined;
    Float64* inCombined;
    fftw_plan pBackward ;
    fftw_complex* precomputedFFT_spcFluxOverErr2;
    fftw_complex* precomputedFFT_spcOneOverErr2;

};


}

#endif

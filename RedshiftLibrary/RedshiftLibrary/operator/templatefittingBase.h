#ifndef _REDSHIFT_OPERATOR_TEMPLATE_FITTING_BASE_
#define _REDSHIFT_OPERATOR_TEMPLATE_FITTING_BASE_

#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/statistics/priorhelper.h>
#include <RedshiftLibrary/spectrum/template/template.h>

#include <vector>

namespace NSEpic
{

class CSpectrum;
class COperatorResult;
class CModelSpectrumResult;

/**
 * \ingroup Redshift
 */
class COperatorTemplateFittingBase : public COperator
{

public:
    //Rule of 5 defaults
    COperatorTemplateFittingBase() = default;
    COperatorTemplateFittingBase(COperatorTemplateFittingBase const& other) = default;
    COperatorTemplateFittingBase(COperatorTemplateFittingBase&& other) = default;
    COperatorTemplateFittingBase& operator=(COperatorTemplateFittingBase const& other) = default;
    COperatorTemplateFittingBase& operator=(COperatorTemplateFittingBase&& other) = default;
    virtual ~COperatorTemplateFittingBase()=default;

    virtual  std::shared_ptr<COperatorResult> Compute( const CSpectrum& spectrum,
                                                       const CTemplate& tpl,
                                                      const TFloat64Range& lambdaRange,
                                                       const TFloat64List& redshifts,
                                                       Float64 overlapThreshold,
                                                       std::vector<CMask> additional_spcMasks,
                                                       std::string opt_interp,
                                                       Int32 opt_extinction,
                                                       Int32 opt_dustFitting,
                                                       CPriorHelper::TPriorZEList logprior=CPriorHelper::TPriorZEList(),
                                                       Bool keepigmism = false,
                                                       Float64 FitEbmvCoeff=-1,
                                                       Float64 FitMeiksinIdx=-1) = 0;

  Int32  ComputeSpectrumModel(const CSpectrum& spectrum,
                              const CTemplate& tpl,
                              Float64 redshift,
                              Float64 EbmvCoeff,
                              Int32 meiksinIdx,
                              Float64 amplitude,
                              std::string opt_interp,
                              std::string opt_extinction,
                              const TFloat64Range& lambdaRange,
                              Float64 overlapThreshold,
                              std::shared_ptr<CModelSpectrumResult> & spc);

  inline virtual bool IsFFTProcessing() {return false;}; 

protected:
  Int32  RebinTemplate( const CSpectrum& spectrum,
                          const CTemplate& tpl, 
                          Float64 redshift,
                          const TFloat64Range& lambdaRange,
                          std::string opt_interp,
                          TFloat64Range& currentRange,
                          Float64& overlaprate,
                          Float64 overlapThreshold);// const;

  CTemplate       m_templateRebined_bf; //buffer
  CSpectrumSpectralAxis m_spcSpectralAxis_restframe; //buffer
  CMask           m_mskRebined_bf; //buffer
  //Likelihood
  Float64 EstimateLikelihoodCstLog(const CSpectrum& spectrum, const TFloat64Range& lambdaRange);

};


}

#endif

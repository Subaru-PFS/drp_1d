#include <RedshiftLibrary/processflow/inputcontext.h>
#include <RedshiftLibrary/processflow/parameterstore.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/template/template.h>

using namespace NSEpic;

CInputContext::CInputContext(std::shared_ptr<CSpectrum> spc,
                             std::shared_ptr<CTemplateCatalog> tmplCatalog,
                             std::shared_ptr<CRayCatalog> rayCatalog,
                             std::shared_ptr<CParameterStore> paramStore):
  m_Spectrum(std::move(spc)),
  m_TemplateCatalog(std::move(tmplCatalog)),
  m_RayCatalog(std::move(rayCatalog)),
  m_ParameterStore(std::move(paramStore))
{
  InitSpectrum();
  std::string calibrationPath =  m_ParameterStore->Get<std::string>( "calibrationDir");  
  InitIsmIgm(calibrationPath);

  TFloat64Range lambdaRange = m_ParameterStore->Get<TFloat64Range>("lambdarange");
  m_Spectrum->GetSpectralAxis().ClampLambdaRange( lambdaRange, m_lambdaRange );
}


void CInputContext::InitSpectrum()
{
    // Smooth flux
    Int64 smoothWidth;
    m_ParameterStore->Get( "smoothWidth", smoothWidth, 0 );
    if( smoothWidth > 0 )
        m_Spectrum->GetFluxAxis().ApplyMeanSmooth(smoothWidth);

    // Continuum removal params
    std::string medianRemovalMethod;
    m_ParameterStore->Get( "continuumRemoval.method", medianRemovalMethod, "IrregularSamplingMedian" );
    m_Spectrum->SetContinuumEstimationMethod(medianRemovalMethod);

    Float64 medianKernelWidth;
    m_ParameterStore->Get( "continuumRemoval.medianKernelWidth", medianKernelWidth, 75.0 );
    m_Spectrum->SetMedianWinsize(medianKernelWidth);

    Float64 nscales;
    m_ParameterStore->Get( "continuumRemoval.decompScales", nscales, 6.0 );
    m_Spectrum->SetDecompScales((Int32)nscales);

    std::string dfBinPath;
    m_ParameterStore->Get( "continuumRemoval.binPath", dfBinPath, "absolute_path_to_df_binaries_here" );
    m_Spectrum->SetWaveletsDFBinPath(dfBinPath);
}

void CInputContext::InitIsmIgm(const std::string & calibrationPath)
{
    
    Float64 ebmv_start=0.0;
    Float64 ebmv_step=0.1;
    UInt32 ebmv_n=10;
    m_ParameterStore->Get( "ebmv.start", ebmv_start, 0. );
    m_ParameterStore->Get( "ebmv.step", ebmv_step, 0.1 );
    m_ParameterStore->Get( "ebmv.count", ebmv_n, 10 );
    //ISM
    auto ismCorrectionCalzetti = std::make_shared<CSpectrumFluxCorrectionCalzetti>();
    ismCorrectionCalzetti->Init(calibrationPath, ebmv_start, ebmv_step, ebmv_n);
    //IGM
    auto igmCorrectionMeiksin = std::make_shared<CSpectrumFluxCorrectionMeiksin>();
    igmCorrectionMeiksin->Init(calibrationPath);

    //push in all templates
    for(std::string s : m_TemplateCatalog->GetCategoryList()){ 
        TTemplateRefList  TplList = m_TemplateCatalog->GetTemplate(TStringList{s});
        for (auto tpl : TplList)
        {
            tpl->m_ismCorrectionCalzetti = ismCorrectionCalzetti;
            if(s!="star")//no igm for stars
                tpl->m_igmCorrectionMeiksin = igmCorrectionMeiksin;
        }   
    }
}


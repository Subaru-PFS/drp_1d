#include <RedshiftLibrary/linemodel/modelspectrumresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <RedshiftLibrary/spectrum/spectrum.h>

using namespace NSEpic;

/**
 * \brief Empty constructor.
 **/
CModelSpectrumResult::CModelSpectrumResult()
{
  
}

/**
 * \brief Sets the model to CSpectrum ( spc ).
 **/
CModelSpectrumResult::CModelSpectrumResult(const CSpectrum& spc):
    m_model(spc)
{
    //probably can add model params as class variable here..
}
CModelSpectrumResult::CModelSpectrumResult(const CSpectrum& spc, 
                                            const std::string tplName,
                                            const Float64 dustCoeff,
                                            const Int32 meiksinIdx,
                                            const Float64 amplitude):
  m_model(spc),
  m_dustCoeff(dustCoeff),
  m_meiksinIdx(meiksinIdx),
  m_tplName(tplName),
  m_amplitude(amplitude)
{
}
/**
 * \brief Empty destructor.
 **/
CModelSpectrumResult::~CModelSpectrumResult()
{

}

/**
 * \brief Prints to argument stream each lambda and flux in the model.
 **/
void CModelSpectrumResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    const CSpectrumSpectralAxis& spectralAxis = m_model.GetSpectralAxis();
    const CSpectrumFluxAxis& modelFluxAxis = m_model.GetFluxAxis();

    stream <<  "#lambda\tflux\t"<< std::endl;
    for ( int i=0; i<spectralAxis.GetSamplesCount(); i++)
    {
        stream << spectralAxis[i] << std::setprecision(16) << "\t" << std::scientific << modelFluxAxis[i] << std::endl;
    }


}

/**
 * \brief Empty method.
 **/
void CModelSpectrumResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}

CSpectrum& CModelSpectrumResult::GetSpectrum()
{
    return m_model;
}

void CModelSpectrumResult::getData(const std::string& name, double **data, int *size) const
{
  if( name.compare("ModelLambda") == 0)
    {
      *size = m_model.GetSpectralAxis().GetSamplesCount();
      *data = const_cast<double *>(m_model.GetSpectralAxis().GetSamples());
    }
  else if( name.compare("ModelFlux") == 0)
    {
      *size = m_model.GetFluxAxis().GetSamplesCount();
      *data = const_cast<double *>(m_model.GetFluxAxis().GetSamples());
    }
  
}

void CModelSpectrumResult::getData(const std::string& name, Int32& v) const
{
  if(name.compare("ModelMeiksinIdx") == 0)
    v = m_meiksinIdx;  
}
void CModelSpectrumResult::getData(const std::string& name, std::string& v) const
{
  if(name.compare("ModelTplName") == 0)
    v = m_tplName;
}
void CModelSpectrumResult::getData(const std::string& name, Float64& v) const
{
  if(name.compare("ModelDustCoeff") == 0)
    v = m_dustCoeff;
  else 
    if(name.compare("ModelAmplitude") == 0)
      v = m_amplitude; 
}

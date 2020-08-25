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
  if( name.compare("model_lambda") == 0)
    {
      *size = model.GetSpectralAxis().GetSamplesCount();
      *data = const_cast<double *>(model.GetSpectralAxis().GetSamples());
    }
  else if( name.compare("model_flux") == 0)
    {
      *size = model.GetFluxAxis().GetSamplesCount();
      *data = const_cast<double *>(model.GetFluxAxis().GetSamples());
    }
  
}

#include <RedshiftLibrary/operator/modelspectrumresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <RedshiftLibrary/spectrum/spectrum.h>

using namespace NSEpic;

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

#include <RedshiftLibrary/operator/modelspectrumresult.h>

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
/*
CModelSpectrumResult::CModelSpectrumResult()
{
  this->m_type = "CModelSpectrumResult";
  ModelLambda = TFloat64List(0);
  ModelFlux = TFloat64List(0);
}
*/
/**
 * \brief Sets the model to CSpectrum ( spc ).
 **/
CModelSpectrumResult::CModelSpectrumResult(const CSpectrum& spc):
  m_model(spc),
  ModelLambda(m_model.GetSpectralAxis().GetSamplesVector()),
  ModelFlux(m_model.GetFluxAxis().GetSamplesVector())
{
  this->m_type = "CModelSpectrumResult";

  
  //probably can add model params as class variable here..
}

/**
 * \brief Empty destructor.
 **/
CModelSpectrumResult::~CModelSpectrumResult()
{

}


CSpectrum& CModelSpectrumResult::GetSpectrum()
{
    return m_model;
}


#include "RedshiftLibrary/operator/modelspectrumresult.h"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include "RedshiftLibrary/spectrum/spectrum.h"

using namespace NSEpic;

/**
 * \brief Sets the model to CSpectrum ( spc ).
 **/
CModelSpectrumResult::CModelSpectrumResult(const CSpectrum& spc):
  ModelLambda(spc.GetSpectralAxis().GetSamplesVector()),
  ModelFlux(spc.GetFluxAxis().GetSamplesVector())
{
  this->m_type = "CModelSpectrumResult";
  //probably can add model params as class variable here..
}

CModelSpectrumResult::CModelSpectrumResult(CSpectrum&& spc):
  ModelLambda(std::move(spc.GetSpectralAxis().GetSamplesVector())),
  ModelFlux(std::move(spc.GetFluxAxis().GetSamplesVector()))
{
  this->m_type = "CModelSpectrumResult";  
};

/*
CModelSpectrumResult::CModelSpectrumResult():
  ModelLambda(TFloat64List(0)),
  ModelFlux(TFloat64List(0))
{}
*/


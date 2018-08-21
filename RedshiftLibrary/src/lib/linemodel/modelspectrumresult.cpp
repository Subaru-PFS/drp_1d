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
CModelSpectrumResult::CModelSpectrumResult(CSpectrum spc)
{
    model = CSpectrum(spc);
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
    const CSpectrumSpectralAxis& spectralAxis = model.GetSpectralAxis();
    const CSpectrumFluxAxis& modelFluxAxis = model.GetFluxAxis();

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
    return model;
}

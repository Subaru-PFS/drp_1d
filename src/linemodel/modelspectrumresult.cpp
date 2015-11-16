#include <epic/redshift/linemodel/modelspectrumresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

#include <epic/redshift/spectrum/spectrum.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CModelSpectrumResult )

CModelSpectrumResult::CModelSpectrumResult()
{
}

CModelSpectrumResult::CModelSpectrumResult(CSpectrum spc)
{
    model = CSpectrum(spc);
}

CModelSpectrumResult::~CModelSpectrumResult()
{

}

Void CModelSpectrumResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    const CSpectrumSpectralAxis& spectralAxis = model.GetSpectralAxis();
    const CSpectrumFluxAxis& modelFluxAxis = model.GetFluxAxis();

    stream <<  "#lambda\tflux\t"<< std::endl;
    for ( int i=0; i<spectralAxis.GetSamplesCount(); i++)
    {
        stream <<  spectralAxis[i] << std::setprecision(16) << "\t" << std::scientific << modelFluxAxis[i] << std::endl;
    }


}

Void CModelSpectrumResult::SaveLine(const CDataStore &store, std::ostream& stream ) const
{

}

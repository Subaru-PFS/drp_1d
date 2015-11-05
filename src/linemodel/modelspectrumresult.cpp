#include <epic/redshift/linemodel/modelspectrumresult.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <fstream>
#include <iomanip>      // std::setprecision

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CModelSpectrumResult )

CModelSpectrumResult::CModelSpectrumResult()
{

}

CModelSpectrumResult::~CModelSpectrumResult()
{

}

Void CModelSpectrumResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
//    CSpectrumSpectralAxis& spectralAxis = model.GetSpectralAxis();
//    CSpectrumFluxAxis& modelFluxAxis = model.GetFluxAxis();

//    stream <<  "#lambda\tflux\t"<< std::endl;
//    for ( int i=0; i<spectralAxis.size(); i++)
//    {
//        stream <<  spectralAxis[i] << std::setprecision(16) << "\t" << std::scientific << modelFluxAxis[i] << std::endl;
//    }


}

Void CModelSpectrumResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{

}

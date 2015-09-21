#ifndef ELEMENT_H
#define ELEMENT_H


#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/ray/catalog.h>

namespace NSEpic
{

class CLineModelElement
{

public:

    CLineModelElement();
    ~CLineModelElement();


    virtual void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift) =0;
    virtual void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )=0;

    virtual Float64 GetFittedAmplitude(Int32 subeIdx)=0;

    Int32 FindElementIndex(Int32 LineCatalogIndex);

protected:

    bool m_OutsideLambdaRange;
    std::vector<Int32> m_LineCatalogIndexes;

private:



};

}







#endif // ELEMENT_H


#ifndef LINEMODEL_ELEMENT_SINGLELINE_H
#define LINEMODEL_ELEMENT_SINGLELINE_H

#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace NSEpic
{

class CSingleLine:public CLineModelElement
{

public:

    CSingleLine(const CRay &r, Float64 nominalWidth, std::vector<Int32> catalogIndexes);
    ~CSingleLine();

    void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift);
    void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    Float64 GetFittedAmplitude(Int32 subeIdx);

private:

    void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift);

    CRay    m_Ray;
    Float64 m_NominalWidth;
    Float64 m_FittedAmplitude;

    Float64 m_NSigmaSupport;
    Int32 m_Start;
    Int32 m_End;

};

}



#endif // ELEMENT_H


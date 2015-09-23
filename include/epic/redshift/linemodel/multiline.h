#ifndef LINEMODEL_ELEMENT_MULTILINE_H
#define LINEMODEL_ELEMENT_MULTILINE_H

#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace NSEpic
{

class CMultiLine:public CLineModelElement
{

public:

    CMultiLine(std::vector<CRay> rs, std::vector<Float64> nominalAmplitudes, Float64 nominalWidth, std::vector<Int32> catalogIndexes);
    ~CMultiLine();


    void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift);
    void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    Float64 GetFittedAmplitude(Int32 subeIdx);
    void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);

private:
    Int32 FindElementIndex(std::string LineTagStr);
    void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift);

    std::vector<CRay>       m_Rays;
    Float64                 m_NominalWidth;
    std::vector<Float64>        m_FittedAmplitudes;
    std::vector<Float64>        m_NominalAmplitudes;

    Float64                     m_NSigmaSupport;
    std::vector<Int32>          m_Start;
    std::vector<Int32>          m_End;

    std::vector<bool>           m_OutsideLambdaRangeList;
};

}


#endif // ELEMENT_H


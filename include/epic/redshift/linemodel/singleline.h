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
    //Float64 FitAmplitudeIterative( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end); //deprecated
    void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    Float64 GetFittedAmplitude(Int32 subeIdx);
    void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);
    bool IsOutsideLambdaRange(Int32 subeIdx);

private:
    Float64 GetLineWidth(Float64 lambda, Float64 z);
    Int32 FindElementIndex(std::string LineTagStr);
    void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift);

    CRay    m_Ray;
    Float64 m_SignFactor;
    Float64 m_NominalWidth;
    Float64 m_FittedAmplitude;

    Float64 m_NSigmaSupport;
    Int32 m_Start;
    Int32 m_End;

};

}



#endif // ELEMENT_H


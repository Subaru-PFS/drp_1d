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

    CSingleLine(const CRay &r, Int32 widthType, Float64 nominalWidth, std::vector<Int32> catalogIndexes);
    ~CSingleLine();

    std::string GetRayName(Int32 subeIdx);

    void prepareSupport(const CSpectrumSpectralAxis& spectralAxis, Float64 redshift);
    TInt32RangeList getSupport();

    void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift);
    //Float64 FitAmplitudeIterative( const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64 lambda, Float64 width, Int32 start, Int32 end); //deprecated
    Float64 getModelAtLambda( Float64 lambda, Float64 redshift );
    void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift );
    void initSpectrumModel( CSpectrumFluxAxis& modelfluxAxis );
    Float64 GetNominalAmplitude(Int32 subeIdx);
    Float64 GetFittedAmplitude(Int32 subeIdx);
    Float64 GetFittedAmplitudeErrorSigma(Int32 subeIdx);
    Float64 GetElementAmplitude();
    void SetFittedAmplitude(Float64 A);
    void LimitFittedAmplitude(Int32 subeIdx, Float64 limit);
    bool IsOutsideLambdaRange(Int32 subeIdx);



private:
    Float64 GetLineWidth(Float64 lambda, Float64 z);
    Int32 FindElementIndex(std::string LineTagStr);

    CRay    m_Ray;
    Float64 m_SignFactor;
    Float64 m_NominalWidth;
    Float64 m_FittedAmplitude;
    Float64 m_FittedAmplitudeErrorSigma;

    Float64 m_NSigmaSupport;
    Int32 m_Start;
    Int32 m_End;

};

}



#endif // ELEMENT_H


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
    enum ELineWidthType
    {
        nWidthType_PSFInstrumentDriven = 1,
        nWidthType_ZDriven = 2,
    };

    CLineModelElement();
    ~CLineModelElement();


    virtual void fitAmplitude(const CSpectrumSpectralAxis& spectralAxis, const CSpectrumFluxAxis& fluxAxis, Float64  redshift) =0;
    virtual void addToSpectrumModel( const CSpectrumSpectralAxis& modelspectralAxis, CSpectrumFluxAxis& modelfluxAxis, Float64 redshift )=0;

    virtual Float64 GetFittedAmplitude(Int32 subeIdx)=0;
    virtual void LimitFittedAmplitude(Int32 subeIdx, Float64 limit)=0;

    Int32 GetSize();
    bool IsOutsideLambdaRange();
    virtual bool IsOutsideLambdaRange(Int32 subeIdx)=0;
    Int32 FindElementIndex(Int32 LineCatalogIndex);
    virtual Int32 FindElementIndex(std::string LineTagStr)=0;

protected:

    Int32 m_LineWidthType;
    Float64 m_Resolution;
    Float64 m_FWHM_factor;

    Float64 m_OutsideLambdaRangeOverlapThreshold;
    bool m_OutsideLambdaRange;
    std::vector<Int32> m_LineCatalogIndexes;

private:



};

}







#endif // ELEMENT_H


#include <epic/redshift/linemodel/element.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>


using namespace NSEpic;

CLineModelElement::CLineModelElement(const std::string& widthType, const Float64 resolution, const Float64 velocity)
{
    m_LineWidthType = widthType;

    //m_Resolution = 250.0 * (1.0 + 0.0); //dr=+0.5 found empirically on VVDS DEEP 651
    m_Resolution = resolution;
    m_Velocity = velocity;
    m_FWHM_factor = 2.35;


    m_NSigmaSupport = 8.0;

    m_OutsideLambdaRange = true;
    m_OutsideLambdaRangeOverlapThreshold = 0.1;
    //example: 0.1 means 10% of the line is allowed to be outside the spectrum with the line still considered inside the lambda range
}

CLineModelElement::~CLineModelElement()
{
}

std::string CLineModelElement::GetElementTypeTag()
{
    return m_ElementType;
}

Int32 CLineModelElement::FindElementIndex(Int32 LineCatalogIndex)
{
    Int32 idx = -1;
    for( UInt32 iElts=0; iElts<m_LineCatalogIndexes.size(); iElts++ )
    {
        if(m_LineCatalogIndexes[iElts] == LineCatalogIndex){
            idx = iElts;
            break;
        }
    }

    return idx;
}

Int32 CLineModelElement::GetSize()
{
    return (Int32)m_LineCatalogIndexes.size();
}

bool CLineModelElement::IsOutsideLambdaRange()
{
    return m_OutsideLambdaRange;
}

Float64 CLineModelElement::GetLineWidth(Float64 lambda, Float64 z, Bool isEmission)
{
    Float64 instrumentSigma = 0.0;
    Float64 velocitySigma = 0.0;
    if( m_LineWidthType == "psfinstrumentdriven"){
        instrumentSigma = lambda/m_Resolution/m_FWHM_factor;
    }else if( m_LineWidthType == "zdriven"){
        instrumentSigma = m_NominalWidth*(1+z);
    }else if( m_LineWidthType == "fixed"){
        instrumentSigma = m_NominalWidth;
    }else if( m_LineWidthType == "fixedvelocity"){
        Float64 v = m_Velocity;
        Float64 c = 299792.458;
        velocitySigma = v/c*lambda*(1+z);//, useless /(1+z)*(1+z);
        if(!isEmission){
            velocitySigma*=3.0;
        }
    }

    Float64 sigma = sqrt(instrumentSigma*instrumentSigma + velocitySigma*velocitySigma);


    return sigma;
}

Float64 CLineModelElement::GetNSigmaSupport()
{
    return m_NSigmaSupport;
}

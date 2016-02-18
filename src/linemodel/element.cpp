#include <epic/redshift/linemodel/element.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>


using namespace NSEpic;

CLineModelElement::CLineModelElement(const std::string& widthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption)
{
    m_LineWidthType = widthType;

    //m_Resolution = 250.0 * (1.0 + 0.0); //dr=+0.5 found empirically on VVDS DEEP 651
    m_Resolution = resolution;
    m_VelocityEmission = velocityEmission;
    m_VelocityAbsorption= velocityAbsorption;
    m_FWHM_factor = 2.35;

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

Float64 CLineModelElement::GetLineWidth(Float64 redshiftedlambda, Float64 z, Bool isEmission)
{
    Float64 instrumentSigma = 0.0;
    Float64 velocitySigma = 0.0;
    if( m_LineWidthType == "psfinstrumentdriven"){
        instrumentSigma = redshiftedlambda/m_Resolution/m_FWHM_factor;
    }else if( m_LineWidthType == "zdriven"){
        instrumentSigma = m_NominalWidth*(1+z);
    }else if( m_LineWidthType == "fixed"){
        instrumentSigma = m_NominalWidth;
    }else if( m_LineWidthType == "fixedvelocity"){
        Float64 v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        Float64 c = 300000.0;
        Float64 pfsSimuCompensationFactor = 1.0;
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;//, useless /(1+z)*(1+z);

        instrumentSigma = redshiftedlambda/m_Resolution/m_FWHM_factor;
    }

    Float64 sigma = sqrt(instrumentSigma*instrumentSigma + velocitySigma*velocitySigma);


    return sigma;
}

Float64 CLineModelElement::GetLineProfile(std::string profile, Float64 xc, Float64 sigma)
{
    Float64 val=0.0;

    if(profile=="SYM"){
        const Float64 xsurc = xc/sigma;
        val = exp(-0.5*xsurc*xsurc);
    }else if(profile=="ASYM"){
        const Float64 xsurc = xc/sigma/2.5;
        const Float64 alpha=10.0;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
    }
    return val;
}

Float64 CLineModelElement::GetLineProfileDerivSigma(std::string profile, Float64 x, Float64 x0, Float64 sigma)
{
    Float64 val=0.0;
    Float64 cel = 300000.0;
    Float64 xc = x-x0;
    if(true || profile=="SYM"){
        const Float64 xsurc = xc/sigma;
        val = xc*xc /cel *x0 /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);
    }else if(profile=="ASYM"){
        //not implemeted yet... to be calculated...
        const Float64 xsurc = xc/sigma/2.5;
        const Float64 alpha=10.0;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
    }
    return val;
}

Float64 CLineModelElement::GetNSigmaSupport(std::string profile)
{
    static Float64 nominal = 8.0;
    Float64 val=nominal;

    if(profile=="SYM"){
        val = nominal;
    }else if(profile=="ASYM"){
        val = nominal*2.5;
    }
    return val;
}

void CLineModelElement::SetVelocityEmission(Float64 vel)
{
    m_VelocityEmission = vel;
}

void CLineModelElement::SetVelocityAbsorption(Float64 vel)
{
    m_VelocityAbsorption = vel;
}

Float64 CLineModelElement::GetVelocityEmission()
{
    return m_VelocityEmission;
}

Float64 CLineModelElement::GetVelocityAbsorption()
{
    return m_VelocityAbsorption;
}

#include <epic/redshift/linemodel/element.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/spectrum.h>
#include <epic/core/log/log.h>

#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <vector>

namespace bfs = boost::filesystem;


using namespace NSEpic;

CLineModelElement::CLineModelElement(const std::string& widthType, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption)
{
    m_LineWidthType = widthType;

    //m_Resolution = 250.0 * (1.0 + 0.0); //dr=+0.5 found empirically on VVDS DEEP 651
    m_Resolution = resolution;
    m_VelocityEmission = velocityEmission;
    m_VelocityAbsorption= velocityAbsorption;
    m_instrumentResolutionEmpiricalFactor = 230.0/325.0/2.35; //derived from (emission line) linemodel-width fit on VUDS ECDFS flags3+4

    m_asym_sigma_coeff = 1.0;
    m_asym_alpha = 4.5;

    m_symxl_sigma_coeff = 5.0;

    m_asym2_sigma_coeff = 2.0;
    m_asym2_alpha = 2.0;

    m_asymfit_sigma_coeff = 2.0;
    m_asymfit_alpha = 2.0;
    m_asymfit_delta = 0.0;

    m_OutsideLambdaRange = true;
    m_OutsideLambdaRangeOverlapThreshold = 0.33; //33% overlap minimum in order to keep the line
    //example: 0.33 means 66% of the line is allowed to be outside the spectrum with the line still considered inside the lambda range

    //LoadDataExtinction(); //uncomment if this line profile is used
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

Float64 CLineModelElement::GetLineWidth(Float64 redshiftedlambda, Float64 z, Bool isEmission, std::string profile)
{
    if(profile=="EXTINCT")
    {
        return 300*(1.0+z); //hardcoded, in Angstrom
    }

    Float64 instrumentSigma = 0.0;
    Float64 velocitySigma = 0.0;

    if( m_LineWidthType == "instrumentdriven"){
        instrumentSigma = redshiftedlambda/m_Resolution*m_instrumentResolutionEmpiricalFactor;
    }else if( m_LineWidthType == "fixed"){
        instrumentSigma = m_NominalWidth;
    }else if( m_LineWidthType == "combined"){
        Float64 v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        Float64 c = 300000.0;
        Float64 pfsSimuCompensationFactor = 1.0;
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;//, useless /(1+z)*(1+z);
        instrumentSigma = redshiftedlambda/m_Resolution*m_instrumentResolutionEmpiricalFactor;
    }else if( m_LineWidthType == "velocitydriven"){
        Float64 v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        Float64 c = 300000.0;
        Float64 pfsSimuCompensationFactor = 1.0;
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;//, useless /(1+z)*(1+z);
    }else if( m_LineWidthType == "nispsim2016"){
        instrumentSigma = (redshiftedlambda*8.121e-4 + 7.4248)/2.35;
        Float64 v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        Float64 c = 300000.0;
        Float64 pfsSimuCompensationFactor = 1.0;
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;
    }

    Float64 sigma = sqrt(instrumentSigma*instrumentSigma + velocitySigma*velocitySigma);


    return sigma;
}

Float64 CLineModelElement::GetLineProfile(std::string profile, Float64 x, Float64 x0, Float64 sigma)
{
    Float64 xc = x-x0;
    Float64 val=0.0;

    if(profile=="SYM"){
        const Float64 xsurc = xc/sigma;
        val = exp(-0.5*xsurc*xsurc);
    }else if(profile=="SYMXL"){
        const Float64 coeff = m_symxl_sigma_coeff;
        sigma = sigma*coeff;
        const Float64 xsurc = xc/sigma;
        val = exp(-0.5*xsurc*xsurc);
    }else if(profile=="LAU"){
        const Float64 xsurc = xc/sigma;
        const Float64 x = xsurc*2.0;
        val = 1.0/(1+x*x);
    }else if(profile=="ASYM"){
        const Float64 coeff = m_asym_sigma_coeff;

        sigma = sigma*coeff;
        const Float64 xsurc = xc/sigma;
        const Float64 alpha = m_asym_alpha;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
    }else if(profile=="ASYM2"){
        const Float64 coeff = m_asym2_sigma_coeff;

        sigma = sigma*coeff;
        const Float64 xsurc = xc/sigma;
        const Float64 alpha = m_asym2_alpha;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
    }else if(profile=="ASYMFIT" || profile.find("ASYMFIXED")!=std::string::npos){
        const Float64 coeff = m_asymfit_sigma_coeff;
        const Float64 xcd = xc+m_asymfit_delta;

        sigma = sigma*coeff;
        const Float64 xsurc = xcd/sigma;
        const Float64 alpha = m_asymfit_alpha;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
    }else if(profile=="EXTINCT"){
        Float64 sigma_rest = m_dataN*m_dataStepLambda;
        Float64 z = sigma/sigma_rest - 1.0;
        Float64 dataStartLambda = (x0/(1+z)) - sigma_rest/2.0;
        Int32 valI = int( (x/(1.0+z)-dataStartLambda)/m_dataStepLambda );
        return m_dataExtinctionFlux[valI];
    }

    return val;
}

Float64 CLineModelElement::GetLineProfileDerivSigma(std::string profile, Float64 x, Float64 x0, Float64 sigma)
{
    Float64 val=0.0;
    Float64 cel = 300000.0;
    Float64 xc = x-x0;
    if(profile=="SYM"){
        const Float64 xsurc = xc/sigma;
        val = xc*xc /cel *x0 /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);
    }else if(profile=="SYMXL"){
        const Float64 coeff = m_symxl_sigma_coeff;
        sigma = sigma*coeff;
        const Float64 xsurc = xc/sigma;
        val = xc*xc /cel *x0 /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);
    }else if(profile=="ASYM"){
        const Float64 coeff = m_asym_sigma_coeff;

        sigma = sigma*coeff;
        const Float64 xsurc = xc/sigma;
        const Float64 alpha = m_asym_alpha;
        const Float64 valsym = exp(-0.5*xsurc*xsurc);
        const Float64 valsymd = xc*xc /cel *x0 /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);

        const Float64 valasym = (1.0+erf(alpha/sqrt(2.0)*xsurc));
        const Float64 arg = alpha*xc/sqrt(2)/sigma;
        //const Float64 valasymd = -alpha/sqrt(2*M_PI)*xc /(sigma*sigma) /cel*x0*exp(-arg*arg);
        const Float64 valasymd = -alpha*sqrt(2)/sqrt(M_PI)*xc /(sigma*sigma) /cel*x0*exp(-arg*arg);
        val = valsym*valasymd+valsymd*valasym;
        //val = valsymd;

        //Float64 v = sigma*cel/x0;
        //val = -sqrt(2)*alpha*cel*(x - x0)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))*exp(-pow(alpha*cel*(x - x0), 2)/(2*pow(v*x0, 2)))/(sqrt(M_PI)*pow(v,2)*x0) + 1.0*pow(cel*(x - x0),2)*(erf(sqrt(2)*alpha*cel*(x - x0)/(2*v*x0)) + 1)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))/(pow(v,3)*pow(x0,2));

    }else if(profile=="ASYM2"){
        const Float64 coeff = m_asym2_sigma_coeff;

        sigma = sigma*coeff;
        const Float64 xsurc = xc/sigma;
        const Float64 alpha = m_asym2_alpha;
        const Float64 valsym = exp(-0.5*xsurc*xsurc);
        const Float64 valsymd = xc*xc /cel *x0 /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);

        const Float64 valasym = (1.0+erf(alpha/sqrt(2.0)*xsurc));
        const Float64 arg = alpha*xc/sqrt(2)/sigma;
        //const Float64 valasymd = -alpha/sqrt(2*M_PI)*xc /(sigma*sigma) /cel*x0*exp(-arg*arg);
        const Float64 valasymd = -alpha*sqrt(2)/sqrt(M_PI)*xc /(sigma*sigma) /cel*x0*exp(-arg*arg);
        val = valsym*valasymd+valsymd*valasym;
        //val = valsymd;

        //Float64 v = sigma*cel/x0;
        //val = -sqrt(2)*alpha*cel*(x - x0)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))*exp(-pow(alpha*cel*(x - x0), 2)/(2*pow(v*x0, 2)))/(sqrt(M_PI)*pow(v,2)*x0) + 1.0*pow(cel*(x - x0),2)*(erf(sqrt(2)*alpha*cel*(x - x0)/(2*v*x0)) + 1)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))/(pow(v,3)*pow(x0,2));

    }else if(profile=="ASYMFIT"  || profile.find("ASYMFIXED")!=std::string::npos){
        const Float64 coeff = m_asymfit_sigma_coeff;
        const Float64 xcd = xc+m_asymfit_delta;

        sigma = sigma*coeff;
        const Float64 xsurc = xcd/sigma;
        const Float64 alpha = m_asymfit_alpha;
        const Float64 valsym = exp(-0.5*xsurc*xsurc);
        const Float64 valsymd = xcd*xcd /cel *x0 /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);

        const Float64 valasym = (1.0+erf(alpha/sqrt(2.0)*xsurc));
        const Float64 arg = alpha*xcd/sqrt(2)/sigma;
        //const Float64 valasymd = -alpha/sqrt(2*M_PI)*xcd /(sigma*sigma) /cel*x0*exp(-arg*arg);
        const Float64 valasymd = -alpha*sqrt(2)/sqrt(M_PI)*xcd /(sigma*sigma) /cel*x0*exp(-arg*arg);
        val = valsym*valasymd+valsymd*valasym;
        //val = valsymd;

        //Float64 v = sigma*cel/x0;
        //val = -sqrt(2)*alpha*cel*(x - x0)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))*exp(-pow(alpha*cel*(x - x0), 2)/(2*pow(v*x0, 2)))/(sqrt(M_PI)*pow(v,2)*x0) + 1.0*pow(cel*(x - x0),2)*(erf(sqrt(2)*alpha*cel*(x - x0)/(2*v*x0)) + 1)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))/(pow(v,3)*pow(x0,2));

    }else if(profile=="EXTINCT"){ //NOT IMPLEMENTED FOR THIS PROFILE, to be done when necessary...
        Float64 sigma_rest = m_dataN*m_dataStepLambda;
        Float64 z = sigma/sigma_rest - 1.0;
        Float64 dataStartLambda = (x0/(1+z)) - sigma_rest/2.0;
        Int32 valI = int( (x/(1.0+z)-dataStartLambda)/m_dataStepLambda );
        return m_dataExtinctionFlux[valI];
    }
    return val;
}

Float64 CLineModelElement::GetNSigmaSupport(std::string profile)
{
    static Float64 nominal = 8.0;
    Float64 val=nominal;

    if(profile=="SYM"){
        val = nominal;
    }else if(profile=="LAU"){
        val = nominal*2.0;
    }else if(profile=="ASYM"){
        val = nominal*m_asym_sigma_coeff;
    }else if(profile=="ASYM2"){
        val = nominal*m_asym2_sigma_coeff;
    }else if(profile=="SYMXL"){
        val = nominal*m_symxl_sigma_coeff;
    }else if(profile=="ASYMFIT"  || profile.find("ASYMFIXED")!=std::string::npos){
        val = nominal*m_asymfit_sigma_coeff*2.5;
    }else if(profile=="EXTINCT"){
        val = 1.0;
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

void CLineModelElement::SetAsymfitWidthCoeff(Float64 coeff)
{
    m_asymfit_sigma_coeff = coeff;
}

Float64 CLineModelElement::GetAsymfitWidthCoeff()
{
    return m_asymfit_sigma_coeff;
}

void CLineModelElement::SetAsymfitAlphaCoeff(Float64 coeff)
{
    m_asymfit_alpha = coeff;
}

Float64 CLineModelElement::GetAsymfitAlphaCoeff()
{
    return m_asymfit_alpha;
}

void CLineModelElement::SetAsymfitDelta(Float64 coeff)
{
    m_asymfit_delta = coeff;
}

Float64 CLineModelElement::GetAsymfitDelta()
{
    return m_asymfit_delta;
}

Float64 CLineModelElement::GetSumCross()
{
    return m_sumCross;
}

Float64 CLineModelElement::GetSumGauss()
{
    return m_sumGauss;
}



/**
 * Laod the extinction residue data.
 */
Bool CLineModelElement::LoadDataExtinction()
{
    std::string filePathStr = "/home/aschmitt/Documents/amazed/methods/linemodel/extinction_element/extinction_data/element_meiksin_dl0.1cropped930-1230.fits.txt";
    Log.LogDebug ( "Parsing ASCII file %s.", filePathStr.c_str() );
    if( !bfs::exists( filePathStr.c_str() ) )
    {
        Log.LogError( "Read: Path for element_extinction file does not exist." );
        return false;
    }

    bfs::ifstream file;
    file.open( filePathStr.c_str() );

    m_dataExtinctionFlux = (Float64 *) malloc(m_dataN* sizeof(Float64));

    Int32 i = 0;
    file.clear();
    file.seekg( 0 );
    int l = file.tellg();
    for( std::string line; std::getline( file, line ); )
    {
        if( !boost::starts_with( line, "#" ) )
        {
            std::istringstream iss( line );
            Float64 x, y;
            iss >> x >> y;
            //spcSpectralAxis[i] = x; //not reading x, hardcoded from 930 to 1230, with dlambda = 0.1A
            m_dataExtinctionFlux[i] = y;
            i++;
            if(i>=m_dataN)
            {
                break;
            }
        }
    }
    file.close();
    if(i<m_dataN)
    {
        Log.LogError( "Read: linemodel- extinction residue: data not read successfully" );
    }
    Log.LogInfo( "Read exctinction successfully" );
    return true;
}

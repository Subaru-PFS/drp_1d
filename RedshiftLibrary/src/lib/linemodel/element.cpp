#include <RedshiftLibrary/linemodel/element.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/log/log.h>

#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <vector>

namespace bfs = boost::filesystem;


using namespace NSEpic;

CLineModelElement::CLineModelElement(const std::string& widthType, const Float64 nsigmasupport, const Float64 resolution, const Float64 velocityEmission, const Float64 velocityAbsorption)
{
    if( widthType == "instrumentdriven"){
        m_LineWidthType = INSTRUMENTDRIVEN;
    } else if( widthType == "fixed"){
        m_LineWidthType = FIXED;
    } else if( widthType == "combined"){
        m_LineWidthType = COMBINED;
    } else if( widthType == "velocitydriven"){
        m_LineWidthType = VELOCITYDRIVEN;
    } else if( widthType == "nispsim2016"){
        m_LineWidthType = NISPSIM2016;
    } else if( widthType == "nispvsspsf201707"){
        m_LineWidthType = NISPVSSPSF201707;
    } else {
        Log.LogError("Unknown LineWidthType %s", widthType.c_str());
        throw std::runtime_error("Unknown LineWidthType");
    }

    m_nsigmasupport = nsigmasupport;

    //m_Resolution = 250.0 * (1.0 + 0.0); //dr=+0.5 found empirically on VVDS DEEP 651
    m_Resolution = resolution;
    m_VelocityEmission = velocityEmission;
    m_VelocityAbsorption= velocityAbsorption;
    m_instrumentResolutionEmpiricalFactor = 230.0/325.0/2.35; //derived from (emission line) linemodel-width fit on VUDS ECDFS flags3+4

    m_SourceSizeDispersion = 0.1;

    m_asym_sigma_coeff = 1.0;
    m_asym_alpha = 4.5;

    m_symxl_sigma_coeff = 5.0;

    m_asym2_sigma_coeff = 2.0;
    m_asym2_alpha = 2.0;

    m_asymfit_sigma_coeff = 2.0;
    m_asymfit_alpha = 2.0;
    m_asymfit_delta = 0.0;

    m_dataExtinctionFlux = NULL;

    m_OutsideLambdaRange = true;
    m_OutsideLambdaRangeOverlapThreshold = 0.33; //33% overlap minimum in order to keep the line
    //example: 0.33 means 66% of the line is allowed to be outside the spectrum with the line still considered inside the lambda range

    //LoadDataExtinction(); //uncomment if this line profile is used
    m_fittingGroupInfo = "-1";
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

Float64 CLineModelElement::GetLineWidth(Float64 redshiftedlambda, Float64 z, Bool isEmission, CRay::TProfile profile)
{
    if(profile==CRay::EXTINCT)
    {
        return 300*(1.0+z); //hardcoded, in Angstrom
    }

    Float64 instrumentSigma = 0.0;
    Float64 velocitySigma = 0.0;
    Float64 sourcesizeSigma = 0.0;

    const Float64 c = 300000.0;
    const Float64 pfsSimuCompensationFactor = 1.0;
    const Float64 arcsecPix = 0.355; //0.3 is the same as in tipsfast
    const Float64 angstromPix = 13.4;
    Float64 v;

    switch (m_LineWidthType) {
    case INSTRUMENTDRIVEN:
        instrumentSigma = redshiftedlambda/m_Resolution*m_instrumentResolutionEmpiricalFactor;
        break;
    case FIXED:
        instrumentSigma = m_NominalWidth;
        break;
    case COMBINED:
        v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;//, useless /(1+z)*(1+z);
        instrumentSigma = redshiftedlambda/m_Resolution*m_instrumentResolutionEmpiricalFactor;
        break;
    case VELOCITYDRIVEN:
        v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;//, useless /(1+z)*(1+z);
        break;
    case NISPSIM2016:
        instrumentSigma = (redshiftedlambda*8.121e-4 + 7.4248)/2.35;
        v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;
        break;
    case NISPVSSPSF201707:
        //+ considers Instrument PSF=f_linearregression(lambda) from MDB-EE50: SpaceSegment.PLM.PLMAsRequired.PLMNISPrEE50rEE80
        //      arcsec/pixel from : SpaceSegment.Instrument.NISP.NISPAsRequired.NISPGRAPSFRefEE50 : (0.355)
        //      angstrom/pixel from : SpaceSegment.Instrument.NISP.NISPAsRequired.NISPGRAAverageDlambda
        //      Leads to linear regression: sigma_psf = 3.939e-4*wl_angstrom + 2.191
        //+ considers source size in the dispersion direction
        //+ considers velocity
        instrumentSigma = (redshiftedlambda*3.939e-4 + 2.191); //probably a realistic calib.
        //instrumentSigma = (redshiftedlambda*4.661e-4 + 2.593); //2017b calib
        //instrumentSigma = 11.; //(approx. 10 or 11?) for tips-fast current version 201708

        sourcesizeSigma = m_SourceSizeDispersion*angstromPix/arcsecPix;

        v = m_VelocityEmission;
        if(!isEmission){
            v = m_VelocityAbsorption;
        }
        velocitySigma = pfsSimuCompensationFactor*v/c*redshiftedlambda;
        break;
    default:
        Log.LogError("Invalid LineWidthType %d", m_LineWidthType);
        throw std::runtime_error("Unknown LineWidthType");
    }

    Float64 sigma = sqrt(instrumentSigma*instrumentSigma + velocitySigma*velocitySigma + sourcesizeSigma*sourcesizeSigma);


    return sigma;
}

Float64 CLineModelElement::GetLineProfile(CRay::TProfile profile, Float64 x, Float64 x0, Float64 sigma)
{
    Float64 xc = x-x0;
    Float64 val = 0.0;
    Float64 xsurc;
    Float64 coeff;
    Float64 alpha;

    Float64 sigma_rest, z, dataStartLambda, xcd;
    Float64 muz, sigmaz, delta, gamma1, m0;
    Int32 valI;

    switch (profile) {
    case CRay::SYM:
        xsurc = xc/sigma;
        val = exp(-0.5*xsurc*xsurc);
        break;
    case CRay::SYMXL:
        coeff = m_symxl_sigma_coeff;
        sigma = sigma*coeff;
        xsurc = xc/sigma;
        val = exp(-0.5*xsurc*xsurc);
        break;
    case CRay::LOR:
        xsurc = xc/sigma;
        val = 1.0/(1+xsurc*xsurc);
        //height of the peak is 2*A/pi/c
        break;
    case CRay::ASYM:
        coeff = m_asym_sigma_coeff;

        sigma = sigma*coeff;
        xsurc = xc/sigma;
        alpha = m_asym_alpha;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
        break;
    case CRay::ASYM2:
        coeff = m_asym2_sigma_coeff;

        sigma = sigma*coeff;
        xsurc = xc/sigma;
        alpha = m_asym2_alpha;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
        break;
    case CRay::ASYMFIT:
    case CRay::ASYMFIXED:
        coeff = m_asymfit_sigma_coeff;
        sigma = sigma*coeff;

        //*
        //correction in order to have the line shifted on the mean: from https://en.wikipedia.org/wiki/Skew_normal_distribution
        delta = m_asymfit_alpha/std::sqrt(1.+m_asymfit_alpha*m_asymfit_alpha);
        muz = delta*sqrt(2./M_PI);
        xc = xc + m_asymfit_sigma_coeff*muz;
        //*/

        /*
        //correction in order to have the line shifted on the mode: from https://en.wikipedia.org/wiki/Skew_normal_distribution
        delta = m_asymfit_alpha/std::sqrt(1.+m_asymfit_alpha*m_asymfit_alpha);
        muz = delta*sqrt(2./M_PI);
        sigmaz = std::sqrt(1-muz*muz);
        gamma1 = ((4-M_PI)/2.0)*pow(delta*std::sqrt(2/M_PI), 3.)/pow(1-2*delta*delta/M_PI, 3./2.);
        m0 = muz - gamma1*sigmaz/2.0 - 0.5*exp(-2*M_PI/m_asymfit_alpha);
        xc = xc + m_asymfit_sigma_coeff*m0;
        //*/

        xcd = xc+m_asymfit_delta;
        xsurc = xcd/sigma;
        alpha = m_asymfit_alpha;
        val = exp(-0.5*xsurc*xsurc)*(1.0+erf(alpha/sqrt(2.0)*xsurc));
        break;
    case CRay::EXTINCT:
        sigma_rest = m_dataN*m_dataStepLambda;
        z = sigma/sigma_rest - 1.0;
        dataStartLambda = (x0/(1+z)) - sigma_rest/2.0;
        valI = int( (x/(1.0+z)-dataStartLambda)/m_dataStepLambda );
        return m_dataExtinctionFlux[valI];
    default:
        Log.LogError("Invalid ray profile %d", profile);
        throw std::runtime_error("Invalid profile");
    }

    //WARNING/TODO/CHECK: this allows multirollmodel to fit the fluxes directly
    //use sigma normalized profiles
    //val /= sigma;

    return val;
}


Float64 CLineModelElement::GetLineFlux(CRay::TProfile profile, Float64 sigma, Float64 A)
{
    Float64 val=0.0;
    switch (profile) {
    case CRay::SYM:
        val = A*sigma*sqrt(2*M_PI);
        break;
    case CRay::LOR:
        val = A*sigma*M_PI;
        break;
    case CRay::ASYM:
        val = A*sigma*sqrt(2*M_PI);
        break;
    case CRay::ASYMFIT:
    case CRay::ASYMFIXED:
        val = A*sigma*m_asymfit_sigma_coeff*sqrt(2*M_PI); //not checked if this analytic integral is correct
        break;
    default:
        Log.LogError("Invalid ray profile for GetLineFlux : %d", profile);
        throw std::runtime_error("Invalid profile");
    }
    return val;
}

Float64 CLineModelElement::GetLineProfileDerivZ(CRay::TProfile profile, Float64 x, Float64 lambda0, Float64 redshift, Float64 sigma){
  Float64 xc = x-lambda0*(1+redshift);
  Float64 val=0.0;
  Float64 xsurc, coeff, alpha, xcd;

  switch (profile) {
  case CRay::SYM:
      xsurc = xc/sigma;
      val = lambda0 /sigma * xsurc * exp(-0.5*xsurc*xsurc);
      break;
  case CRay::SYMXL:
      coeff = m_symxl_sigma_coeff;
      sigma = sigma*coeff;
      xsurc = xc/sigma;
      val = lambda0 /sigma * xsurc * exp(-0.5*xsurc*xsurc);
      break;
  case CRay::ASYM:
      coeff = m_asym_sigma_coeff;
      sigma = sigma*coeff;
      xsurc = xc/sigma;
      alpha = m_asym_alpha;
      val = lambda0 /sigma * xsurc * exp(-0.5*xsurc*xsurc) *(1.0+erf(alpha/sqrt(2.0)*xsurc)) -alpha * lambda0 /sqrt(2*M_PI) /sigma * exp(-(1+alpha*alpha)/2 * xsurc * xsurc);
      break;
  case CRay::ASYM2:
      coeff = m_asym2_sigma_coeff;
      sigma = sigma*coeff;
      xsurc = xc/sigma;
      alpha = m_asym2_alpha;
      val = lambda0 /sigma * xsurc * exp(-0.5*xsurc*xsurc) *(1.0+erf(alpha/sqrt(2.0)*xsurc)) -alpha * lambda0 /sqrt(2*M_PI) /sigma * exp(-(1+alpha*alpha)/2 * xsurc * xsurc);
      break;
  case CRay::ASYMFIT:
  case CRay::ASYMFIXED:
      coeff = m_asymfit_sigma_coeff;
      xcd = xc+m_asymfit_delta;

      sigma = sigma*coeff;
      xsurc = xcd/sigma;
      alpha = m_asymfit_alpha;
      val = lambda0 /sigma * xsurc * exp(-0.5*xsurc*xsurc) *(1.0+erf(alpha/sqrt(2.0)*xsurc)) -alpha * lambda0 /sqrt(2*M_PI) /sigma * exp(-(1+alpha*alpha)/2 * xsurc * xsurc);
      break;
  default:
      Log.LogError("Deriv for Z not IMPLEMENTED for profile %d", profile);
      throw std::runtime_error("Invalid profile");
  }
  return val;
}

Float64 CLineModelElement::GetLineProfileDerivVel(CRay::TProfile profile, Float64 x, Float64 x0, Float64 sigma, Bool isEmission){
    const Float64 c = 300000.0;
    const Float64 pfsSimuCompensationFactor = 1.0;
    Float64 v, v_to_sigma;

    switch (m_LineWidthType) {
    case INSTRUMENTDRIVEN:
    case FIXED:
        return 0.0;
    case COMBINED:
    case NISPSIM2016:
    case NISPVSSPSF201707: //not supported as of 2017-07
        v = isEmission ? m_VelocityEmission : m_VelocityAbsorption;
        v_to_sigma = pfsSimuCompensationFactor/c*x0; //velocity sigma = v_to_sigma * v
        return v_to_sigma * v_to_sigma * v /sigma * GetLineProfileDerivSigma(profile, x, x0, sigma);
    case VELOCITYDRIVEN:
        v_to_sigma = pfsSimuCompensationFactor/c*x0;
        return v_to_sigma * GetLineProfileDerivSigma(profile, x, x0, sigma);
    default:
        Log.LogError("Invalid LineWidthType : %d", m_LineWidthType);
        throw std::runtime_error("Unknown LineWidthType");
    }
    return 0.0;
}

Float64 CLineModelElement::GetLineProfileDerivSigma(CRay::TProfile profile, Float64 x, Float64 x0, Float64 sigma)
{
    Float64 val=0.0;
    //Float64 cel = 300000.0;
    Float64 xc = x-x0;
    Float64 xsurc, coeff, alpha, valsym, valsymd, valasym, arg, valasymd, xcd;
    Float64 sigma_rest, z, dataStartLambda;
    Int32 valI;

    switch (profile) {
    case CRay::SYM:
        xsurc = xc/sigma;
        val = xc*xc  /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);
        break;
    case CRay::SYMXL:
        coeff = m_symxl_sigma_coeff;
        sigma = sigma*coeff;
        xsurc = xc/sigma;
        val = xc*xc/(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);
        break;
    case CRay::ASYM:
        coeff = m_asym_sigma_coeff;

        sigma = sigma*coeff;
        xsurc = xc/sigma;
        alpha = m_asym_alpha;
        valsym = exp(-0.5*xsurc*xsurc);
        valsymd = xc*xc/(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);

        valasym = (1.0+erf(alpha/sqrt(2.0)*xsurc));
        arg = alpha*xc/sqrt(2)/sigma;
        //const Float64 valasymd = -alpha/sqrt(2*M_PI)*xc /(sigma*sigma) /cel*x0*exp(-arg*arg);
        valasymd = -alpha*sqrt(2)/sqrt(M_PI)*xc /(sigma*sigma)*exp(-arg*arg);
        val = valsym*valasymd+valsymd*valasym;
        //val = valsymd;

        //Float64 v = sigma*cel/x0;
        //val = -sqrt(2)*alpha*cel*(x - x0)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))*exp(-pow(alpha*cel*(x - x0), 2)/(2*pow(v*x0, 2)))/(sqrt(M_PI)*pow(v,2)*x0) + 1.0*pow(cel*(x - x0),2)*(erf(sqrt(2)*alpha*cel*(x - x0)/(2*v*x0)) + 1)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))/(pow(v,3)*pow(x0,2));
        break;
    case CRay::ASYM2:
        coeff = m_asym2_sigma_coeff;

        sigma = sigma*coeff;
        xsurc = xc/sigma;
        alpha = m_asym2_alpha;
        valsym = exp(-0.5*xsurc*xsurc);
        valsymd = xc*xc /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);

        valasym = (1.0+erf(alpha/sqrt(2.0)*xsurc));
        arg = alpha*xc/sqrt(2)/sigma;
        //const Float64 valasymd = -alpha/sqrt(2*M_PI)*xc /(sigma*sigma) /cel*x0*exp(-arg*arg);
        valasymd = -alpha*sqrt(2)/sqrt(M_PI)*xc /(sigma*sigma)*exp(-arg*arg);
        val = valsym*valasymd+valsymd*valasym;
        //val = valsymd;

        //Float64 v = sigma*cel/x0;
        //val = -sqrt(2)*alpha*cel*(x - x0)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))*exp(-pow(alpha*cel*(x - x0), 2)/(2*pow(v*x0, 2)))/(sqrt(M_PI)*pow(v,2)*x0) + 1.0*pow(cel*(x - x0),2)*(erf(sqrt(2)*alpha*cel*(x - x0)/(2*v*x0)) + 1)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))/(pow(v,3)*pow(x0,2));
        break;
    case CRay::ASYMFIT:
    case CRay::ASYMFIXED:
        coeff = m_asymfit_sigma_coeff;
        xcd = xc+m_asymfit_delta;

        sigma = sigma*coeff;
        xsurc = xcd/sigma;
        alpha = m_asymfit_alpha;
        valsym = exp(-0.5*xsurc*xsurc);
        valsymd = xcd*xcd  /(sigma*sigma*sigma) * exp(-0.5*xsurc*xsurc);

        valasym = (1.0+erf(alpha/sqrt(2.0)*xsurc));
        arg = alpha*xcd/sqrt(2)/sigma;
        //const Float64 valasymd = -alpha/sqrt(2*M_PI)*xcd /(sigma*sigma) /cel*x0*exp(-arg*arg);
        valasymd = -alpha*sqrt(2)/sqrt(M_PI)*xcd /(sigma*sigma)*exp(-arg*arg);
        val = valsym*valasymd+valsymd*valasym;
        //val = valsymd;

        //Float64 v = sigma*cel/x0;
        //val = -sqrt(2)*alpha*cel*(x - x0)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))*exp(-pow(alpha*cel*(x - x0), 2)/(2*pow(v*x0, 2)))/(sqrt(M_PI)*pow(v,2)*x0) + 1.0*pow(cel*(x - x0),2)*(erf(sqrt(2)*alpha*cel*(x - x0)/(2*v*x0)) + 1)*exp(-0.5*pow(cel*(x - x0)/(v*x0),2))/(pow(v,3)*pow(x0,2));
        break;
    case CRay::EXTINCT:
        //NOT IMPLEMENTED FOR THIS PROFILE, to be done when necessary...
        Float64 sigma_rest = m_dataN*m_dataStepLambda;
        Float64 z = sigma/sigma_rest - 1.0;
        Float64 dataStartLambda = (x0/(1+z)) - sigma_rest/2.0;
        Int32 valI = int( (x/(1.0+z)-dataStartLambda)/m_dataStepLambda );
        return m_dataExtinctionFlux[valI];
    }
    return val;
}

Float64 CLineModelElement::GetNSigmaSupport(CRay::TProfile profile)
{
    Float64 nominal = m_nsigmasupport;
    Float64 val=nominal;

    switch (profile) {
    case CRay::SYM:
        val = nominal;
        break;
    case CRay::LOR:
        val = nominal*2.0;
        break;
    case CRay::ASYM:
        val = nominal*m_asym_sigma_coeff;
        break;
    case CRay::ASYM2:
        val = nominal*m_asym2_sigma_coeff;
        break;
    case CRay::SYMXL:
        val = nominal*m_symxl_sigma_coeff;
        break;
    case CRay::ASYMFIT:
    case CRay::ASYMFIXED:
        val = nominal*m_asymfit_sigma_coeff*2.5;
        break;
    case CRay::EXTINCT:
        val = 1.0;
        break;
    default:
        Log.LogError("Invalid ray profile for GetNSigmaSupport : %d", profile);
        throw std::runtime_error("Invalid profile");
    }
    return val;
}


void CLineModelElement::SetSourcesizeDispersion(Float64 sigma)
{
    m_SourceSizeDispersion = sigma;
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

Float64 CLineModelElement::GetVelocity()
{
    Float64 vel=-1;
    if(m_Rays.size()>0)
    {
        if(m_Rays[0].GetIsEmission())
        {
            vel = m_VelocityEmission;
        }else{
            vel = m_VelocityAbsorption;
        }
    }
    return vel;
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

void CLineModelElement::SetSumCross(Float64 val)
{
    m_sumCross=val;
}

Float64 CLineModelElement::GetDtmFree()
{
    return m_dtmFree;
}

void CLineModelElement::SetDtmFree(Float64 val)
{
    m_dtmFree=val;
}

Float64 CLineModelElement::GetSumGauss()
{
    return m_sumGauss;
}

void CLineModelElement::SetSumGauss(Float64 val)
{
    m_sumGauss=val;
}

Float64 CLineModelElement::GetFitAmplitude()
{
    return m_fitAmplitude;
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

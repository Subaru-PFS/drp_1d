#include "RedshiftLibrary/spectrum/LSFbyExtrinsicComponents.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;
using namespace std;


CLSFbyExtrinsicComponents::~CLSFbyExtrinsicComponents()
{}

/**
 * Note that LSFType correspond to the instrument responses (including the source response)
 * Note that FROMSPECTRUMDATA includes this latter
*/
CLSFbyExtrinsicComponents::CLSFbyExtrinsicComponents(const std::string LSFType, 
                                                     const Float64 resolution, 
                                                     const Float64 nominalWidth):
m_Resolution(resolution),
m_NominalWidth(nominalWidth)
{
    //"FIXED"/"NISPSIM2016"/"NISPVSSPSF201707/FROMSPECTRUMDATA
        
    if( LSFType == "fromspectrumdata"){
        m_type = FROMSPECTRUMDATA;
    } else if( LSFType == "fixed"){
        m_type = FIXED;
    } else if( LSFType == "nispsim2016"){
        m_type = NISPSIM2016;
    } else if( LSFType == "nispvsspsf201707"){
        m_type = NISPVSSPSF201707;
    } else {
        Log.LogError("Unknown LSFType %s", LSFType.c_str());
        throw std::runtime_error("Unknown LSFType");
    }
    m_SourceSizeDispersion = 0.1;
    m_instrumentResolutionEmpiricalFactor = 230.0/325.0/2.35; //derived from (emission line) linemodel-width fit on VUDS ECDFS flags3+4

}
void CLSFbyExtrinsicComponents::SetSourcesizeDispersion(Float64 sigma) 
{
    m_SourceSizeDispersion = sigma;
}

//using extrinsic params only
Float64 CLSFbyExtrinsicComponents::GetWidth(Float64 lambda) const
{
    const Float64 arcsecPix = 0.355; //0.3 is the same as in tipsfast
    const Float64 angstromPix = 13.4;
    Float64 instrumentSigma, sourcesizeSigma = 0.0;        
    Float64 defaultSigma = lambda/m_Resolution*m_instrumentResolutionEmpiricalFactor;
    switch(m_type)
    {
        case FROMSPECTRUMDATA: //this corresponds to reading info from LSF, assuming that an LSF is defined
            
            instrumentSigma = (m_width==0.)?defaultSigma:m_width;//could be replaced by a call to a function
            break;
        case FIXED:
            instrumentSigma = m_NominalWidth;//which should be read from param.json
            break;
        case NISPSIM2016:
            instrumentSigma = (lambda*8.121e-4 + 7.4248)/2.35;
            break;
        case NISPVSSPSF201707:
            //+ considers Instrument PSF=f_linearregression(lambda) from MDB-EE50: SpaceSegment.PLM.PLMAsRequired.PLMNISPrEE50rEE80
            //      arcsec/pixel from : SpaceSegment.Instrument.NISP.NISPAsRequired.NISPGRAPSFRefEE50 : (0.355)
            //      angstrom/pixel from : SpaceSegment.Instrument.NISP.NISPAsRequired.NISPGRAAverageDlambda
            //      Leads to linear regression: sigma_psf = 3.939e-4*wl_angstrom + 2.191
            //+ considers source size in the dispersion direction
            //+ considers velocity
            instrumentSigma = (lambda*3.939e-4 + 2.191); //probably a realistic calib.
            //instrumentSigma = (lambda*4.661e-4 + 2.593); //2017b calib
            //instrumentSigma = 11.; //(approx. 10 or 11?) for tips-fast current version 201708
            sourcesizeSigma = m_SourceSizeDispersion*angstromPix/arcsecPix;
            break;
        default:
            Log.LogError("Invalid LSFType %d", m_type);
            throw std::runtime_error("Unknown LSFType");
    }
    Float64 width = sqrt(instrumentSigma*instrumentSigma + sourcesizeSigma*sourcesizeSigma);
    //m_width = width;
    return width;
}


void CLSFbyExtrinsicComponents::SetWidth(const Float64 width)
{
    m_width = width;
}


bool CLSFbyExtrinsicComponents::IsValid() const
{
    return (m_width>0.0);
}

#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include <RedshiftLibrary/log/log.h>


using namespace NSEpic;
using namespace std;

CLSFGaussianNISPVSSPSF201707::CLSFGaussianNISPVSSPSF201707(Float64 sourcesize):
    CLSF(GaussianNISPVSSPSF201707, std::make_shared<CLineProfileSYM>()),
    m_SourceSizeDispersion(sourcesize)
{
    
}

Float64 CLSFGaussianNISPVSSPSF201707::GetWidth(Float64 lambda) const
{
    const Float64 arcsecPix = 0.355;
    const Float64 angstromPix = 13.4;
    //+ considers Instrument PSF=f_linearregression(lambda) from MDB-EE50: SpaceSegment.PLM.PLMAsRequired.PLMNISPrEE50rEE80
    //      arcsec/pixel from : SpaceSegment.Instrument.NISP.NISPAsRequired.NISPGRAPSFRefEE50 : (0.355)
    //      angstrom/pixel from : SpaceSegment.Instrument.NISP.NISPAsRequired.NISPGRAAverageDlambda
    //      Leads to linear regression: sigma_psf = 3.939e-4*wl_angstrom + 2.191
    //+ considers source size in the dispersion direction
    //+ considers velocity
    Float64 instrumentSigma = (lambda*3.939e-4 + 2.191); //probably a realistic calib.
    //instrumentSigma = (lambda*4.661e-4 + 2.593); //2017b calib
    //instrumentSigma = 11.; //(approx. 10 or 11?) for tips-fast current version 201708
    Float64 sourcesizeSigma = m_SourceSizeDispersion*angstromPix/arcsecPix;

    Float64 width = sqrt(instrumentSigma*instrumentSigma + sourcesizeSigma*sourcesizeSigma);
    return width;
}

bool CLSFGaussianNISPVSSPSF201707::IsValid() const
{
    return true;//TODO
}
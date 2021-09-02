// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/log/log.h"


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
    return true;
}
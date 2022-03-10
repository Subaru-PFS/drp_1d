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
#ifndef _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_
#define _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CSpectrumFluxCorrectionMeiksin
{

public:

    struct MeiksinCorrection{
        TFloat64List lbda; // wavelength
        std::vector<TFloat64List> fluxcorr; // 7 flux correction lists
    };

    CSpectrumFluxCorrectionMeiksin();
    ~CSpectrumFluxCorrectionMeiksin();

    bool LoadCurvesinIncreasingExtinctionOrder( const char* filePath );
    bool Init( std::string calibrationPath, const std::shared_ptr<const CLSF>& lsf, TFloat64Range& lambdaRange);
    TFloat64List Convolve(const TFloat64List& arr, const TFloat64List& kernel);
    TFloat64List ApplyAdaptativeKernel(const TFloat64List& arr, 
                                        const Float64 z_center, 
                                        const std::shared_ptr<const CLSF>& lsf,
                                        const TFloat64List& lambdas);
    void ConvolveAll(const std::shared_ptr<const CLSF>& lsf);
    
    Int32 GetIdxCount() const;
    Int32 GetRedshiftIndex(Float64 z) const;

    TFloat64List GetSegmentsStartRedshiftList() const;
    TFloat64List         GetLSFProfileVector(Float64 lambda0_rest, 
                                             Float64 z_bin_meiksin, 
                                             const std::shared_ptr<const CLSF>& lsf);//for convolution
    Float64 GetLambdaMin() const;
    Float64 GetLambdaMax() const;
    bool meiksinInitFailed = false;
    TFloat64Range m_convolRange;
    std::vector<MeiksinCorrection> m_corrections;
private:
    std::vector<MeiksinCorrection> m_rawCorrections;
    Float64 m_LambdaMin;
    Float64 m_LambdaMax;
    TFloat64List m_kernel;

};

inline TFloat64List CSpectrumFluxCorrectionMeiksin::GetSegmentsStartRedshiftList() const
{
    TFloat64List zstartlist = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0,
                                       4.5, 5.0, 5.5, 6.0, 6.5};
    return zstartlist;
};

inline Float64 CSpectrumFluxCorrectionMeiksin::GetLambdaMin() const
{
    return m_LambdaMin;
}

inline Float64 CSpectrumFluxCorrectionMeiksin::GetLambdaMax() const 
{
    return m_LambdaMax;
}

inline Int32 CSpectrumFluxCorrectionMeiksin::GetIdxCount() const
{
    return 7; //harcoded value from the number of cols in the ascii files
}


}

#endif

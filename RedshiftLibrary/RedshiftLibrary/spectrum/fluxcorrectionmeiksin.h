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

    Bool LoadCurvesinIncreasingExtinctionOrder( const char* filePath );
    Bool Init( std::string calibrationPath, const std::shared_ptr<const CLSF>& lsf);
    TFloat64List Convolve(const TFloat64List& arr, const TFloat64List& kernel);
    TFloat64List ApplyAdaptativeKernel(const TFloat64List& arr, 
                                        const Float64 z_center, 
                                        const std::shared_ptr<const CLSF>& lsf,
                                        const TFloat64List& lambdas);
    void ConvolveAll(const std::shared_ptr<const CLSF>& lsf);
    
    Int32 GetIdxCount() const;
    Int32 GetRedshiftIndex(Float64 z) const;

    std::vector<Float64> GetSegmentsStartRedshiftList() const;
    TFloat64List         GetLSFProfileVector(Float64 lambda0_rest, 
                                             Float64 z_bin_meiksin, 
                                             const std::shared_ptr<const CLSF>& lsf);//for convolution
    Float64 GetLambdaMin() const;
    Float64 GetLambdaMax() const;

    bool meiksinInitFailed = false;

    std::vector<MeiksinCorrection> m_corrections;
private:
    std::vector<MeiksinCorrection> m_rawCorrections;
    Float64 m_LambdaMin;
    Float64 m_LambdaMax;
    TFloat64List m_kernel;

};

inline std::vector<Float64> CSpectrumFluxCorrectionMeiksin::GetSegmentsStartRedshiftList() const
{
    std::vector<Float64> zstartlist = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0,
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

#ifndef _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_
#define _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

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
        std::vector<Float64> lbda; // wavelength
        std::vector<TFloat64List> fluxcorr; // 7 flux correction lists
    };

    CSpectrumFluxCorrectionMeiksin();
    ~CSpectrumFluxCorrectionMeiksin();

    Bool Init( std::string calibrationPath );
    Bool LoadCurvesinIncreasingExtinctionOrder( const char* filePath );

    std::vector<MeiksinCorrection> m_corrections;

    Int32 GetIdxCount() const;
    Int32 GetRedshiftIndex(Float64 z) const;
    std::vector<Float64> GetSegmentsStartRedshiftList() const;

    Float64 GetLambdaMin() const;
    Float64 GetLambdaMax() const;

    bool meiksinInitFailed = false;

private:

    Float64 m_LambdaMin;
    Float64 m_LambdaMax;

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

#ifndef _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_
#define _REDSHIFT_SPECTRUM_FLUXCORRECTIONMEIKSIN_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

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
    Bool LoadFile( const char* filePath );

    std::vector<MeiksinCorrection> m_corrections;

    Int32 GetIdxCount();
    Int32 GetRedshiftIndex(Float64 z);
    std::vector<Float64> GetSegmentsStartRedshiftList();

    Float64  getCoeff(Int32 meiksinIdx, Float64 redshift, Float64 restLambda);


    Float64 GetLambdaMin();
    Float64 GetLambdaMax();

    bool meiksinInitFailed = false;

private:

    Float64 m_LambdaMin;
    Float64 m_LambdaMax;

};

inline std::vector<Float64> CSpectrumFluxCorrectionMeiksin::GetSegmentsStartRedshiftList()
{
    std::vector<Float64> zstartlist = {0.0, 2.0, 2.5, 3.0, 3.5, 4.0,
                                       4.5, 5.0, 5.5, 6.0, 6.5};
    return zstartlist;
};

inline Float64 CSpectrumFluxCorrectionMeiksin::GetLambdaMin()
{
    return m_LambdaMin;
}

inline Float64 CSpectrumFluxCorrectionMeiksin::GetLambdaMax()
{
    return m_LambdaMax;
}

inline Int32 CSpectrumFluxCorrectionMeiksin::GetIdxCount()
{
    return 7; //harcoded value from the number of cols in the ascii files
}


}

#endif

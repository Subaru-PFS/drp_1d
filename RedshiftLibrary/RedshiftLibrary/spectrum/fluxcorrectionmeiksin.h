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
 * /ingroup Redshift

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

    const Float64*  getMeiksinCoeff(Int32 meiksinIdx, Float64 redshift, Float64 maxLambda);


    Float64 GetLambdaMin();
    Float64 GetLambdaMax();

private:

    Float64 m_LambdaMin;
    Float64 m_LambdaMax;

};


}

#endif

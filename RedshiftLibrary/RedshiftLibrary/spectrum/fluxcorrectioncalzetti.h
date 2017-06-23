#ifndef _REDSHIFT_SPECTRUM_FLUXCORRECTIONCALZETTI_
#define _REDSHIFT_SPECTRUM_FLUXCORRECTIONCALZETTI_

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
class CSpectrumFluxCorrectionCalzetti
{

public:
    CSpectrumFluxCorrectionCalzetti();
    ~CSpectrumFluxCorrectionCalzetti();

    Bool Init( std::string calibrationPath );
    Bool LoadFile( const char* filePath );

    Float64 GetLambdaMin();
    Float64 GetLambdaMax();
    Int32 GetNPrecomputedDustCoeffs();

    Float64 GetEbmvValue(Int32 k);

    Float64 getDustCoeff(Int32 kDust, Float64 restLambda );

    const Float64*  getDustCoeff(Float64 dustCoeff, Float64 maxLambda);

    Float64 *m_dataCalzetti;
    Float64 m_NdataCalzetti;

    Int32 m_nDustCoeff;
    Float64 m_dustCoeffStep;
    Float64 m_dustCoeffStart;
    Float64* m_dataDustCoeff;
    bool calzettiInitFailed;

private:

    Float64 m_LambdaMin;
    Float64 m_LambdaMax;



};


}

#endif

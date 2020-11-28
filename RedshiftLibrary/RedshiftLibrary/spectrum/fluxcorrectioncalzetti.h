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
 * \ingroup Redshift
 */
class CSpectrumFluxCorrectionCalzetti
{

public:
    CSpectrumFluxCorrectionCalzetti();
    ~CSpectrumFluxCorrectionCalzetti();

    Bool Init( std::string calibrationPath, Float64 ebmv_start, Float64 ebmv_step, Float64 ebmv_n );
    Bool LoadFile( const char* filePath );

    Float64 GetLambdaMin() const;
    Float64 GetLambdaMax() const;
    Int32 GetNPrecomputedDustCoeffs() const;

    Float64 GetEbmvValue(Int32 k) const;

    Float64 getDustCoeff(Int32 kDust, Float64 restLambda ) const;

    std::vector<Float64> m_dataCalzetti;

    Int32 m_nDustCoeff = 0;
    Float64 m_dustCoeffStep;
    Float64 m_dustCoeffStart;
    std::vector<Float64> m_dataDustCoeff;
    bool calzettiInitFailed = false;

private:

    Float64 m_LambdaMin;
    Float64 m_LambdaMax;



};


}

#endif

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
    Int32 GetNPrecomputedEbmvCoeffs() const;

    Float64 GetEbmvValue(Int32 k) const;

    Float64 GetDustCoeff(Int32 kDust, Float64 restLambda ) const;
    Int32 GetEbmvIndex(Float64 ebmv) const;

    std::vector<Float64> m_dataCalzetti;

    Int32 m_nEbmvCoeff = 0;
    Float64 m_EbmvCoeffStep;
    Float64 m_EbmvCoeffStart;
    std::vector<Float64> m_dataDustCoeff;
    bool calzettiInitFailed = false;

private:

    Float64 m_LambdaMin;
    Float64 m_LambdaMax;



};

}

#endif

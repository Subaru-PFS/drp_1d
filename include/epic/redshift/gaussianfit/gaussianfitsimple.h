#ifndef _REDSHIFT_GAUSSIANFIT_GAUSSIANFIT_
#define _REDSHIFT_GAUSSIANFIT_GAUSSIANFIT_

#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace NSEpic
{

class CSpectrum;

class CGaussianFitSimple
{

public:

    enum EStatus
    {
        nStatus_Success = (1 << 0),

        // Results corresponding to an error
        nStatus_IllegalInput = (1 << 1),
        nStatus_InvalidStudyRange = (1 << 2),
        nStatus_IterationHasNotConverged = (1 << 3),

        // Results corresponding to a success
        nStatus_FailToReachTolerance = (1 << 4)
    };

    CGaussianFitSimple( );
    ~CGaussianFitSimple();

    EStatus    Compute( const CSpectrum& s, const TInt32Range& studyRange );

    Void    GetResults( Float64& amplitude, Float64& position, Float64& width ) const;
    Void    GetResultsPolyCoeff0( Float64& coeff0) const;
    Void    GetResultsError( Float64& amplitude, Float64& position, Float64& width ) const;

private:


    static int GaussF (const gsl_vector *param, void *data, gsl_vector * f);
    static int GaussDF (const gsl_vector * param, void *data, gsl_matrix * J);
    static int GaussFDF (const gsl_vector * param, void *data, gsl_vector * f, gsl_matrix * J);

    Void ComputeFirstGuess( const CSpectrum& spectrum, const TInt32Range& studyRange, Int32 polyOrder, Float64& peakValue, Float64& peakPos, Float64& gaussAmp );

    struct SUserData
    {
        const CSpectrum*    spectrum;
        const TInt32Range*  studyRange;
        Int32               polyOrder;
    };

    Float64 m_Amplitude;
    Float64 m_Mu;
    Float64 m_C;
    Float64 m_AmplitudeErr;
    Float64 m_MuErr;
    Float64 m_CErr;

    Float64 m_coeff0;

    Int32   m_PolyOrder;
    Float64 m_AbsTol;
    Float64 m_RelTol;
};

}

#endif

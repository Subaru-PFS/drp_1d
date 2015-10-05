#ifndef _REDSHIFT_GAUSSIANFIT_MULTIGAUSSIANFIT_
#define _REDSHIFT_GAUSSIANFIT_MULTIGAUSSIANFIT_

#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>


#include <epic/redshift/linemodel/elementlist.h>

#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace NSEpic
{

class CSpectrum;

class CMultiGaussianFit
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

    CMultiGaussianFit( );
    ~CMultiGaussianFit();

    Int32    Compute(CLineModelElementList model );

private:


    static double  my_f (const gsl_vector *v, void *params);

    struct SUserData
    {
        CLineModelElementList* model;
        Int32 nddl;
        std::vector<Int32> modelIdx;
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

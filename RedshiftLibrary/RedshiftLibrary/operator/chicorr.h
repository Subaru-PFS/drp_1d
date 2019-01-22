#ifndef CHICORR_H
#define CHICORR_H

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>
#include <RedshiftLibrary/spectrum/fluxaxis.h>
#include <RedshiftLibrary/spectrum/spectralaxis.h>
#include <RedshiftLibrary/operator/operator.h>


namespace NSEpic
{

class CSpectrum;
class CTemplate;
class CCorrelationResult;
class CChisquareResult;

/**
 * \ingroup Redshift
 */
class COperatorChicorr
{

public:

    COperatorChicorr();
    ~COperatorChicorr();

    int Compute(const CSpectrum& spectrum, const CSpectrum& spectrumWithoutCont, const CTemplate& tpl, const CTemplate& tplWithoutCont, const TFloat64Range& r, const TFloat64List& redhisfts, Float64 overlap, CCorrelationResult *result_corr, CChisquareResult *result_chi);

    Float64 GetComputationDuration() const;

private:

    Float64                 m_TotalDuration;
};


}

#endif // CHICORR_H


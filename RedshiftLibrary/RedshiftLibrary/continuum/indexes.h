#ifndef _REDSHIFT_CONTINUUM_INDEXES_
#define _REDSHIFT_CONTINUUM_INDEXES_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <gsl/gsl_vector.h>
#include <boost/format.hpp>

namespace NSEpic
{
/**
   * \ingroup Redshift
   * \brief
   */
class CContinuumIndexes
{
public:

    struct SContinuumIndex
    {
        Float64 Color;
        Float64 Break;
    };
    struct SContinuumRelevance
    {
        Float64 StdSpectrum;
        Float64 StdContinuum;
    };
    typedef std::vector<SContinuumIndex> TContinuumIndexList;


    TContinuumIndexList getIndexes(const CSpectrum& spectrum, Float64 z);
    SContinuumRelevance getRelevance(const CSpectrum& spectrum, const CSpectrum& continuum);
private:
};

}
#endif

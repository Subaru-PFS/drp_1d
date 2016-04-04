#ifndef _REDSHIFT_CONTINUUM_INDEXES_
#define _REDSHIFT_CONTINUUM_INDEXES_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/spectrum.h>
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
    typedef std::vector<SContinuumIndex> TContinuumIndexList;


    TContinuumIndexList getIndexes(const CSpectrum& spectrum, Float64 z);
private:
    TContinuumIndexList m_indexes;
};

}
#endif

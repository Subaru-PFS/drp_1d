#ifndef _REDSHIFT_CONTINUUM_INDEXES_PRIOR
#define _REDSHIFT_CONTINUUM_INDEXES_PRIOR

#include "RedshiftLibrary/common/datatypes.h"
#include <gsl/gsl_vector.h>
#include <boost/format.hpp>

namespace NSEpic
{
/**
   * \ingroup Redshift
   * \brief
   */
class CContinuumIndexesPrior
{
public:
    CContinuumIndexesPrior();
    bool Init(std::string calibrationPath);
    Float64 GetHeatmapVal( Int32 _index, Float64 _color, Float64 _break);

    typedef std::vector<std::vector<Float64>> TContinuumIndexData;

private:

    Int32 m_nIndexes = 0;

    Float64 m_tbl_color_step = 0;
    Float64 m_tbl_color_min = 0;
    Int32 m_tbl_color_n = 0;

    Float64 m_tbl_break_step = 0;
    Float64 m_tbl_break_min = 0;
    UInt32 m_tbl_break_n = 0;

    std::vector<TContinuumIndexData> m_ciprior_table;
    std::vector<Float64> m_ciprior_max;
};

}
#endif

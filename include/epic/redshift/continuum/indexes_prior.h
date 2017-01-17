#ifndef _REDSHIFT_CONTINUUM_INDEXES_PRIOR
#define _REDSHIFT_CONTINUUM_INDEXES_PRIOR

#include <epic/core/common/datatypes.h>
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

    bool Init(std::string calibrationPath);
    Float64 GetHeatmapVal( Int32 _index, Float64 _color, Float64 _break);

    typedef std::vector<std::vector<Float64>> TContinuumIndexData;

private:

    Int32 m_nIndexes;

    Float64 m_tbl_color_step;
    Float64 m_tbl_color_min;
    Int32 m_tbl_color_n;

    Float64 m_tbl_break_step;
    Float64 m_tbl_break_min;
    Int32 m_tbl_break_n;

    std::vector<TContinuumIndexData> m_ciprior_table;
};

}
#endif

#include "RedshiftLibrary/spectrum/LSFVariableWidth.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/common/exception.h"
using namespace NSEpic;
using namespace std;

//used to create mapping between width and lambda values
//Instead of this 
CLSFGaussianVariableWidth::CLSFGaussianVariableWidth(const std::shared_ptr<const TLSFGaussianVarWidthArgs>& args):
    CLSF(GaussianVariableWidth, std::make_shared<CLineProfileSYM>())
{
    m_width = args->width;
    m_spcAxis = args->lambdas;
    IsValid();
}
Float64 CLSFGaussianVariableWidth::GetWidth(Float64 lambda) const
{
    Int32 idx = -1;
    if(lambda<m_spcAxis[0] || lambda>m_spcAxis[m_width.size()-1])
    {
        Log.LogError("CLSFGaussianVariableWidth::GetWidth: lambda outside spectralAxis range");
        throw std::runtime_error("CLSFGaussianVariableWidth::GetWidth: lambda outside spectralAxis range");
    }

    TFloat64Index::getClosestLowerIndex(m_spcAxis.GetSamplesVector(), lambda, idx);

    if(m_spcAxis[idx]!=lambda)
    {//interpolation
        Float64 t = (lambda - m_spcAxis[idx])/(m_spcAxis[idx+1] - m_spcAxis[idx]);
        Float64 width_interp = m_width[idx] + (m_width[idx+1] - m_width[idx])*t;
        return width_interp;
    }
    return m_width[idx];
}

bool CLSFGaussianVariableWidth::IsValid() const
{
    if(!m_width.size()){
        Log.LogError("CLSFGaussianVariableWidth:: Width array cannot be null");
        throw GlobalException(BAD_COUNTMATCH, "Size is null ");
    }
    if(m_spcAxis.GetSamplesCount() != m_width.size()) {
        Log.LogError("CLSFGaussianVariableWidth:: Spectral axis size and width axis size do not match");
        throw GlobalException(BAD_COUNTMATCH, "Size do not match ");
    };
    for(Float64 w : m_width)
        if(w <= 0. ) return false;
    return true;
}

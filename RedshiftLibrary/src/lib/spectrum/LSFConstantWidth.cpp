#include <RedshiftLibrary/spectrum/LSFConstantWidth.h>
#include <RedshiftLibrary/log/log.h>


using namespace NSEpic;
using namespace std;

/**
 * Constructor.
 */
CLSFGaussianConstantWidth::CLSFGaussianConstantWidth(const Float64 width):
    CLSF(GaussianConstantWidth, std::make_shared<CLineProfileSYM>()),
    m_width(width)
{}

Float64 CLSFGaussianConstantWidth::GetWidth(Float64 lambda) const
{
    return m_width;
}

bool CLSFGaussianConstantWidth::IsValid() const
{
    return (m_width>0.0);
}

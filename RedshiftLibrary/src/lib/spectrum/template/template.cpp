#include <RedshiftLibrary/spectrum/template/template.h>

#include <RedshiftLibrary/common/mask.h>
#include <fstream>

using namespace NSEpic;
using namespace std;

/**
 * Constructor, empty.
 */
CTemplate::CTemplate( )
{

}

/**
 * Constructor, assigns values to members.
 */
CTemplate::CTemplate( const std::string& name, const std::string& category ) :
    m_Category( category ),
    m_Name( name )
{

}

CTemplate::CTemplate( const std::string& name, const std::string& category,
		      CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis) :
    m_Category( category ),
    m_Name( name )
{
  m_SpectralAxis = spectralAxis;
  m_FluxAxis = fluxAxis;
}

/**
 * Destructor, empty.
 */
CTemplate::~CTemplate()
{

}

/**
 * Returns the value stored in m_Name.
 */
const std::string& CTemplate::GetName() const
{
    return m_Name;
}

/**
 * Returns the value stored in m_Category.
 */
const std::string& CTemplate::GetCategory() const
{
    return m_Category;
}

/**
 * Saves the template in the given filePath.
 */
Bool CTemplate::Save( const char* filePath ) const
{
    std::fstream file;

    file.open( filePath, fstream::out );
    if( file.rdstate() & ios_base::failbit )
    {
        return false;
    }

    const CSpectrumSpectralAxis& spectralAxis = GetSpectralAxis();
    CSpectrumSpectralAxis spectralAxisCopy(spectralAxis);
    bool logScale = spectralAxisCopy.IsInLogScale();
    //alway save in Linear Scale
    if(logScale)
    {
        spectralAxisCopy.ConvertToLinearScale();
    }
    const CSpectrumFluxAxis& fluxAxis = GetFluxAxis();
    for ( Int32 i=0; i<GetSampleCount(); i++)
    {
        file.precision(10);
        file  <<  spectralAxisCopy[i] << "\t" ;

        file.precision(10);
        file<< fluxAxis[i] << std::endl;
    }
    file.close();
    return true;
}

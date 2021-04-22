#include <RedshiftLibrary/spectrum/template/template.h>

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <fstream>
#include <iostream>
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
  m_Category( category )
{
    m_Name = name;
}

CTemplate::CTemplate( const std::string& name, const std::string& category,
		      CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis) :
    CSpectrum(spectralAxis, fluxAxis),
    m_Category( category )
{
    m_Name = name;
}

CTemplate::CTemplate( const CTemplate& other): 
    CSpectrum(other),
    m_kDust(other.m_kDust),
    m_meiksinIdx(other.m_meiksinIdx),
    m_Category( other.m_Category),
    m_IsmIgm_kstart(other.m_IsmIgm_kstart),
    m_IsmIgm_kend(other.m_IsmIgm_kend),
    m_computedDustCoeff(other.m_computedDustCoeff), 
    m_computedMeiksingCoeff(other.m_computedDustCoeff),
    m_ismCorrectionCalzetti(other.m_ismCorrectionCalzetti),
    m_igmCorrectionMeiksin(other.m_igmCorrectionMeiksin)
{
        if(other.m_NoIsmIgmFluxAxis.GetSamplesCount())
        m_NoIsmIgmFluxAxis = other.m_NoIsmIgmFluxAxis;
}

CTemplate::CTemplate( const CTemplate& other, TFloat64List mask): 
    CSpectrum(other, mask),
    m_kDust(other.m_kDust),
    m_meiksinIdx(other.m_meiksinIdx),
    m_Category( other.m_Category),
    m_IsmIgm_kstart(other.m_IsmIgm_kstart),
    m_IsmIgm_kend(other.m_IsmIgm_kend),
    m_computedDustCoeff(other.m_computedDustCoeff), 
    m_computedMeiksingCoeff(other.m_computedDustCoeff),
    m_ismCorrectionCalzetti(other.m_ismCorrectionCalzetti),
    m_igmCorrectionMeiksin(other.m_igmCorrectionMeiksin)
{
    if(other.m_NoIsmIgmFluxAxis.GetSamplesCount())
        other.m_NoIsmIgmFluxAxis.MaskAxis(mask, m_NoIsmIgmFluxAxis);
            
}
CTemplate& CTemplate::operator=(const CTemplate& other)
{
    CSpectrum::operator =(other);
    m_Name = other.m_Name;

    if(other.m_NoIsmIgmFluxAxis.GetSamplesCount())
        m_NoIsmIgmFluxAxis = other.m_NoIsmIgmFluxAxis;
        
    m_kDust = other.m_kDust;
    m_meiksinIdx = other.m_meiksinIdx;
    m_computedDustCoeff = other.m_computedDustCoeff; 
    m_computedMeiksingCoeff = other.m_computedDustCoeff;
    m_Category = other.m_Category;
    m_IsmIgm_kstart = other.m_IsmIgm_kstart;
    m_IsmIgm_kend = other.m_IsmIgm_kend;
    m_ismCorrectionCalzetti = other.m_ismCorrectionCalzetti;
    m_igmCorrectionMeiksin = other.m_igmCorrectionMeiksin;
    return *this;
}

/**
 * Returns the value stored in m_Category.
 */
const std::string& CTemplate::GetCategory() const
{
    return m_Category;
}

void CTemplate::SetIsmIgmLambdaRange(TFloat64Range& lbdaRange)
{
    lbdaRange.getClosedIntervalIndices(m_SpectralAxis.GetSamplesVector(), m_IsmIgm_kstart, m_IsmIgm_kend);
}

void CTemplate::SetIsmIgmLambdaRange(Int32 kstart, Int32 kend)
{
    m_IsmIgm_kstart = kstart; m_IsmIgm_kend = kend;
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

//Calzetti extinction
/**
 * if kDust == -1: reset the fluxAxis eliminating only the Dust correction
 */
bool CTemplate::ApplyDustCoeff(Int32 kDust)
{
    if(m_kDust == kDust)
        return true; 
    m_kDust = kDust;

    CSpectrumFluxAxis & FluxAxis = GetFluxAxis();

    for(Int32 k =m_IsmIgm_kstart; k < m_IsmIgm_kend + 1; k++)
    {
        if(m_kDust > -1)
            m_computedDustCoeff[k] = m_ismCorrectionCalzetti->GetDustCoeff( kDust, m_SpectralAxis[k]); 
        else
            m_computedDustCoeff[k] = 1; 
        
        FluxAxis[k] = m_NoIsmIgmFluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
    }
    return true;
}

/**
 * if meiksinIdx == -1: reset the fluxAxis eliminating only the Meiksin correction
 */
bool CTemplate::ApplyMeiksinCoeff(Int32 meiksinIdx, Float64 redshift)
{

    if(m_meiksinIdx == meiksinIdx && m_redshiftMeiksin == redshift)
        return true;
        
    Int32 redshiftIdx = m_igmCorrectionMeiksin->GetRedshiftIndex(redshift); //index for IGM Meiksin redshift range
    m_meiksinIdx = meiksinIdx;
    m_redshiftMeiksin = redshift;
    Bool igmCorrectionAppliedOnce = false;

    CSpectrumFluxAxis & FluxAxis = GetFluxAxis();

    for(Int32 k = m_IsmIgm_kstart; k < m_IsmIgm_kend + 1; k++)
    {
        if(m_SpectralAxis[k] <= m_igmCorrectionMeiksin->GetLambdaMax()){
            if(m_meiksinIdx > -1){
                Int32 kLbdaMeiksin = 0;
                if(m_SpectralAxis[k] >= m_igmCorrectionMeiksin->GetLambdaMin())
                {
                    kLbdaMeiksin = Int32(m_SpectralAxis[k] - m_igmCorrectionMeiksin->GetLambdaMin());
                }else //if lambda lower than min meiksin value, use lower meiksin value
                {
                    kLbdaMeiksin = 0;
                }
                m_computedMeiksingCoeff[k] = m_igmCorrectionMeiksin->m_corrections[redshiftIdx].fluxcorr[meiksinIdx][kLbdaMeiksin];
            }
            else 
                m_computedMeiksingCoeff[k] = 1;
            
            FluxAxis[k] = m_NoIsmIgmFluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
            igmCorrectionAppliedOnce = true;
        }
    }
    return igmCorrectionAppliedOnce;
}

bool CTemplate::CalzettiInitFailed() const
{
    bool failed = false;
    if (!m_ismCorrectionCalzetti)
    {
        failed=true;
    }
    else if (m_ismCorrectionCalzetti->calzettiInitFailed)
    {
        failed = true;
    }
    
    return failed;
}

bool CTemplate::MeiksinInitFailed() const
{
    bool failed = false;
    
    if (!m_igmCorrectionMeiksin)
    {
        failed = true;
    }
    else if (m_igmCorrectionMeiksin->meiksinInitFailed)
    {
        failed = true;
    }

    return failed;
}

//init ism/igm configuration when we change redshift value
bool CTemplate::InitIsmIgmConfig( const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti,
                                const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin)
{
    if(!m_ismCorrectionCalzetti)
        m_ismCorrectionCalzetti = ismCorrectionCalzetti;
    if(!m_igmCorrectionMeiksin)
        m_igmCorrectionMeiksin = igmCorrectionMeiksin;
    return InitIsmIgmConfig();
}

bool CTemplate::InitIsmIgmConfig()
{
    m_kDust = -1;
    m_meiksinIdx = -1;
    m_redshiftMeiksin = -1; 
    m_NoIsmIgmFluxAxis = GetFluxAxis();

    m_computedMeiksingCoeff.resize(m_SpectralAxis.GetSamplesCount());
    std::fill(m_computedMeiksingCoeff.begin(), m_computedMeiksingCoeff.end(), 1.0);
    
    m_computedDustCoeff.resize(m_SpectralAxis.GetSamplesCount());
    std::fill(m_computedDustCoeff.begin(), m_computedDustCoeff.end(), 1.0);
    return true;
}

void CTemplate::ScaleFluxAxis(Float64 amplitude){

    CSpectrum::ScaleFluxAxis(amplitude);
    if (!m_NoIsmIgmFluxAxis.isEmpty())
        m_NoIsmIgmFluxAxis *= amplitude;
}

#include "RedshiftLibrary/spectrum/template/template.h"

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/mask.h"
#include "RedshiftLibrary/spectrum/template/catalog.h"
#include <fstream>
#include <iostream>

#include <numeric>
using namespace NSEpic;
using namespace std;

/**
 * Constructor, assigns values to members.
 */
CTemplate::CTemplate( const std::string& name, const std::string& category ) :
  m_Category( category )
{
    m_Name = name;
}

CTemplate::CTemplate( const std::string& name, const std::string& category,
		      CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis) :
    CSpectrum(std::move(spectralAxis), std::move(fluxAxis)),
    m_Category( category )
{
    m_Name = name;
}

CTemplate::CTemplate( const CTemplate& other): 
    CSpectrum(other),
    m_kDust(other.m_kDust),
    m_meiksinIdx(other.m_meiksinIdx),
    m_meiksinRedshiftIdx(other.m_meiksinRedshiftIdx),
    m_Category( other.m_Category),
    m_IsmIgm_kstart(other.m_IsmIgm_kstart),
    m_Ism_kend(other.m_Ism_kend),
    m_Igm_kend(other.m_Igm_kend),
    m_computedDustCoeff(other.m_computedDustCoeff), 
    m_computedMeiksingCoeff(other.m_computedMeiksingCoeff),
    m_ismCorrectionCalzetti(other.m_ismCorrectionCalzetti),
    m_igmCorrectionMeiksin(other.m_igmCorrectionMeiksin),
    m_NoIsmIgmFluxAxis(other.m_NoIsmIgmFluxAxis)
{
}

CTemplate::CTemplate( CTemplate&& other): 
    CSpectrum(std::move(other)),
    m_kDust(other.m_kDust),
    m_meiksinIdx(other.m_meiksinIdx),
    m_meiksinRedshiftIdx(other.m_meiksinRedshiftIdx),
    m_Category( std::move(other.m_Category)),
    m_IsmIgm_kstart(other.m_IsmIgm_kstart),
    m_Ism_kend(other.m_Ism_kend),
    m_Igm_kend(other.m_Igm_kend),
    m_computedDustCoeff(std::move(other.m_computedDustCoeff)), 
    m_computedMeiksingCoeff(std::move(other.m_computedMeiksingCoeff)),
    m_ismCorrectionCalzetti(std::move(other.m_ismCorrectionCalzetti)),
    m_igmCorrectionMeiksin(std::move(other.m_igmCorrectionMeiksin)),
    m_NoIsmIgmFluxAxis(std::move(other.m_NoIsmIgmFluxAxis))
{
}

CTemplate::CTemplate( const CTemplate& other, const TFloat64List& mask): 
    CSpectrum(other, mask),
    m_kDust(other.m_kDust),
    m_meiksinIdx(other.m_meiksinIdx),
    m_meiksinRedshiftIdx(other.m_meiksinRedshiftIdx),
    m_Category( other.m_Category),
    m_ismCorrectionCalzetti(other.m_ismCorrectionCalzetti),
    m_igmCorrectionMeiksin(other.m_igmCorrectionMeiksin)
{
    if(other.CheckIsmIgmEnabled()){
        TFloat64Range otherRange(other.m_SpectralAxis[other.m_IsmIgm_kstart], other.m_SpectralAxis[other.m_Ism_kend]);
        bool ret = otherRange.getClosedIntervalIndices(m_SpectralAxis.GetSamplesVector(), m_IsmIgm_kstart, m_Ism_kend);
        if(!ret)//complete range is masked
        {
            m_IsmIgm_kstart = -1; m_Ism_kend = -1; m_Igm_kend = -1;//not necessary but for security
            m_NoIsmIgmFluxAxis.clear();
        }else{
            other.m_NoIsmIgmFluxAxis.MaskAxis(mask, m_NoIsmIgmFluxAxis);

            CSpectrumAxis::maskVector(mask, other.m_computedDustCoeff, m_computedDustCoeff);
            CSpectrumAxis::maskVector(mask, other.m_computedMeiksingCoeff, m_computedMeiksingCoeff);
        }
    }
}

CTemplate& CTemplate::operator=(const CTemplate& other)
{
    CSpectrum::operator=(other);

    m_NoIsmIgmFluxAxis = other.m_NoIsmIgmFluxAxis;
    m_kDust = other.m_kDust;
    m_meiksinIdx = other.m_meiksinIdx;
    m_meiksinRedshiftIdx = other.m_meiksinRedshiftIdx;
    m_computedDustCoeff = other.m_computedDustCoeff; 
    m_computedMeiksingCoeff = other.m_computedMeiksingCoeff;
    m_Category = other.m_Category;
    m_IsmIgm_kstart = other.m_IsmIgm_kstart;
    m_Ism_kend = other.m_Ism_kend;
    m_Igm_kend = other.m_Igm_kend;
    m_ismCorrectionCalzetti = other.m_ismCorrectionCalzetti;
    m_igmCorrectionMeiksin = other.m_igmCorrectionMeiksin;
    return *this;
}

CTemplate& CTemplate::operator=(CTemplate&& other)
{
    CSpectrum::operator=(std::move(other));

    m_NoIsmIgmFluxAxis = std::move(other.m_NoIsmIgmFluxAxis);
    m_kDust = other.m_kDust;
    m_meiksinIdx = other.m_meiksinIdx;
    m_meiksinRedshiftIdx = other.m_meiksinRedshiftIdx;
    m_computedDustCoeff = std::move(other.m_computedDustCoeff); 
    m_computedMeiksingCoeff = std::move(other.m_computedMeiksingCoeff);
    m_Category = std::move(other.m_Category);
    m_IsmIgm_kstart = other.m_IsmIgm_kstart;
    m_Ism_kend = other.m_Ism_kend;
    m_Igm_kend = other.m_Igm_kend;
    m_ismCorrectionCalzetti = std::move(other.m_ismCorrectionCalzetti);
    m_igmCorrectionMeiksin = std::move(other.m_igmCorrectionMeiksin);
    return *this;
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

//Calzetti extinction
/**
 * if kDust == -1: reset the fluxAxis eliminating only the Dust correction
 */
bool CTemplate::ApplyDustCoeff(Int32 kDust)
{
    if (!CheckIsmIgmEnabled() || CalzettiInitFailed()){
        Log.LogError("CTemplate::ApplyDustCoeff: try to apply dust extinction without ism initialization");
        throw runtime_error("CTemplate::ApplyDustCoeff: try to apply dust extinction without ism initialization");
    }

    if(m_kDust == kDust)
        return true; 
    m_kDust = kDust;

    CSpectrumFluxAxis & FluxAxis = GetFluxAxis_();

    for(Int32 k =m_IsmIgm_kstart; k < m_Ism_kend + 1; k++)
    {
        if(m_kDust > -1)
            m_computedDustCoeff[k] = m_ismCorrectionCalzetti->GetDustCoeff( kDust, m_SpectralAxis[k]); 
        else
            m_computedDustCoeff[k] = 1.0; 
        
        FluxAxis[k] = m_NoIsmIgmFluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
    }
    return true;
}

/**
 * if meiksinIdx == -1: reset the fluxAxis eliminating only the Meiksin correction
 */
bool CTemplate::ApplyMeiksinCoeff(Int32 meiksinIdx)
{
    if (!CheckIsmIgmEnabled() || MeiksinInitFailed()){
        Log.LogError("CTemplate::ApplyMeiksinCoeff: try to apply igm extinction without igm initialization");
        throw runtime_error("CTemplate::ApplyMeiksinCoeff: try to apply igm extinction without igm initialization");
    }

    if(m_meiksinIdx == meiksinIdx)
        return m_Igm_kend==-1 ? false:true;

    m_meiksinIdx = meiksinIdx;

    if ( m_Igm_kend==-1) 
        return false;
        
    CSpectrumFluxAxis & FluxAxis = GetFluxAxis_();

    for(Int32 k = m_IsmIgm_kstart; k <= m_Igm_kend; k++)
    {
        if(m_meiksinIdx > -1){
            Int32 kLbdaMeiksin = 0;
            if(m_SpectralAxis[k] >= m_igmCorrectionMeiksin->GetLambdaMin())
            {
                kLbdaMeiksin = Int32(m_SpectralAxis[k] - m_igmCorrectionMeiksin->GetLambdaMin());
            }else //if lambda lower than min meiksin value, use lower meiksin value
            {
                kLbdaMeiksin = 0;
            }
            m_computedMeiksingCoeff[k] = m_igmCorrectionMeiksin->m_corrections[m_meiksinRedshiftIdx].fluxcorr[m_meiksinIdx][kLbdaMeiksin];
        }
        else 
            m_computedMeiksingCoeff[k] = 1.0;
        
        FluxAxis[k] = m_NoIsmIgmFluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
    }
    return true;
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
void CTemplate::InitIsmIgmConfig( Float64 redshift,
                                  const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti,
                                  const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin)
{
    InitIsmIgmConfig(0, GetSampleCount()-1, redshift, ismCorrectionCalzetti, igmCorrectionMeiksin);
}

void CTemplate::InitIsmIgmConfig( const TFloat64Range & lbdaRange, Float64 redshift,
                                  const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti,
                                  const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin)
{
    Int32 kstart, kend;
    lbdaRange.getClosedIntervalIndices(m_SpectralAxis.GetSamplesVector(), kstart, kend);
    InitIsmIgmConfig(kstart, kend, redshift, ismCorrectionCalzetti, igmCorrectionMeiksin);
}

void CTemplate::InitIsmIgmConfig( Int32 kstart, Int32 kend, Float64 redshift,
                                  const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti,
                                  const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin)
{
    if (ismCorrectionCalzetti)
        m_ismCorrectionCalzetti = ismCorrectionCalzetti;

    if (igmCorrectionMeiksin)
        m_igmCorrectionMeiksin = igmCorrectionMeiksin;

    if (MeiksinInitFailed() && CalzettiInitFailed() ) {
        Log.LogError("CTemplate::InitIsmIgmConfig: Cannot init ismigm");
        throw runtime_error("CTemplate::InitIsmIgmConfig: Cannot init ismigm");
    }
    
    m_kDust = -1;
    m_meiksinIdx = -1;

    m_IsmIgm_kstart = kstart; 
    m_Ism_kend = kend;
    m_Igm_kend = -1;

    if (!MeiksinInitFailed())
    {
        m_meiksinRedshiftIdx = m_igmCorrectionMeiksin->GetRedshiftIndex(redshift); //index for IGM Meiksin redshift range
    
        // get last index in spectral axis where igm can be applied
        TAxisSampleList::iterator istart = m_SpectralAxis.GetSamplesVector().begin()+m_IsmIgm_kstart;
        TAxisSampleList::iterator iend = m_SpectralAxis.GetSamplesVector().begin()+m_Ism_kend+1;
        TAxisSampleList::iterator it = std::upper_bound(istart, iend, m_igmCorrectionMeiksin->GetLambdaMax());
        if (it!=istart)  // should -1 if not applicable (lambdamax< lmabd[kstart]) 
            m_Igm_kend = it - m_SpectralAxis.GetSamplesVector().begin() -1; 
    }

    if(m_NoIsmIgmFluxAxis.isEmpty()) // initialize when called for the first time
        m_NoIsmIgmFluxAxis = GetFluxAxis();
    else // reset the fluxAxis
        GetFluxAxis_() = m_NoIsmIgmFluxAxis; //note: the type component (raw/continuum/wocontinuum should not have changed)
        
    m_computedMeiksingCoeff.resize(m_SpectralAxis.GetSamplesCount());
    std::fill(m_computedMeiksingCoeff.begin(), m_computedMeiksingCoeff.end(), 1.0);
    
    m_computedDustCoeff.resize(m_SpectralAxis.GetSamplesCount());
    std::fill(m_computedDustCoeff.begin(), m_computedDustCoeff.end(), 1.0);
}

void CTemplate::ScaleFluxAxis(Float64 amplitude){

    CSpectrum::ScaleFluxAxis(amplitude);
    if (CheckIsmIgmEnabled()) 
        m_NoIsmIgmFluxAxis *= amplitude;
}

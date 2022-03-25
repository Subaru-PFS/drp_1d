// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/spectrum/template/template.h"

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/exception.h"
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
        bool ret = otherRange.getClosedIntervalIndices(m_SpectralAxis.GetSamplesVector(), m_IsmIgm_kstart, m_Ism_kend, false);
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
bool CTemplate::Save( const char* filePath ) const
{
    std::fstream file;

    file.open( filePath, fstream::out );
    if( file.rdstate() & ios_base::failbit )
    {
        return false;
    }

    CSpectrumSpectralAxis spectralAxis = GetSpectralAxis();
    bool logScale = spectralAxis.IsInLogScale();
    //alway save in Linear Scale
    if(logScale)
    {
        spectralAxis.ConvertToLinearScale();
    }
    const CSpectrumFluxAxis& fluxAxis = GetFluxAxis();
    for ( Int32 i=0; i<GetSampleCount(); i++)
    {
        file.precision(10);
        file  <<  as_const(spectralAxis)[i] << "\t" ;

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
        throw GlobalException(INTERNAL_ERROR,"CTemplate::ApplyDustCoeff: try to apply dust extinction without ism initialization");
    }

    if(m_kDust == kDust)
        return true; 
    m_kDust = kDust;

    CSpectrumFluxAxis & FluxAxis = GetFluxAxis_();
    const CSpectrumSpectralAxis & SpectralAxis = m_SpectralAxis;
    const CSpectrumFluxAxis & NoIsmIgmFluxAxis = m_NoIsmIgmFluxAxis;

    for(Int32 k =m_IsmIgm_kstart; k < m_Ism_kend + 1; k++)
    {
        if(m_kDust > -1)
            m_computedDustCoeff[k] = m_ismCorrectionCalzetti->GetDustCoeff( kDust, SpectralAxis[k]); 
        else
            m_computedDustCoeff[k] = 1.0; 
        
        FluxAxis[k] = NoIsmIgmFluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
    }
    return true;
}

/**
 * if meiksinIdx == -1: reset the fluxAxis eliminating only the Meiksin correction
 */
bool CTemplate::ApplyMeiksinCoeff(Int32 meiksinIdx)
{
    if (!CheckIsmIgmEnabled() || MeiksinInitFailed()){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::ApplyMeiksinCoeff: try to apply igm extinction without igm initialization");
    }

    if(m_meiksinIdx == meiksinIdx)
        return m_Igm_kend==-1 ? false:true;

    m_meiksinIdx = meiksinIdx;

    if ( m_Igm_kend==-1) 
        return false;
        
    CSpectrumFluxAxis & FluxAxis = GetFluxAxis_();
    const auto & SpectralAxis = m_SpectralAxis;
    const auto & NoIsmIgmFluxAxis = m_NoIsmIgmFluxAxis;

    for(Int32 k = m_IsmIgm_kstart; k <= m_Igm_kend; k++)
    {
        if(m_meiksinIdx > -1){
            Int32 kLbdaMeiksin = 0;
            if(m_SpectralAxis[k] >= m_igmCorrectionMeiksin->GetLambdaMin())
            {
                kLbdaMeiksin = Int32(SpectralAxis[k] - m_igmCorrectionMeiksin->GetLambdaMin());
            }else //if lambda lower than min meiksin value, use lower meiksin value
            {
                kLbdaMeiksin = 0;
            }
            m_computedMeiksingCoeff[k] = m_igmCorrectionMeiksin->m_corrections[m_meiksinRedshiftIdx].fluxcorr[m_meiksinIdx][kLbdaMeiksin];
        }
        else 
            m_computedMeiksingCoeff[k] = 1.0;
        
        FluxAxis[k] = NoIsmIgmFluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
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
    bool ret = lbdaRange.getClosedIntervalIndices(m_SpectralAxis.GetSamplesVector(), kstart, kend);
    if (!ret){
      throw GlobalException(INTERNAL_ERROR,"CTemplate::InitIsmIgmConfig: lambda range outside spectral axis");
    }
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
        throw GlobalException(INTERNAL_ERROR,"CTemplate::InitIsmIgmConfig: Cannot init ismigm");
    }
    
    if (kstart<0 || kstart>=m_SpectralAxis.GetSamplesCount()){
      throw GlobalException(INTERNAL_ERROR,"CTemplate::InitIsmIgmConfig: kstart outside range");
    }
    if (kend<0 || kend>=m_SpectralAxis.GetSamplesCount()){
      throw GlobalException(INTERNAL_ERROR,"CTemplate::InitIsmIgmConfig: kend outside range");
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
        m_Igm_kend = GetIgmEndIndex(m_IsmIgm_kstart, m_Ism_kend);
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

Int32 CTemplate::GetIgmEndIndex(Int32 kstart, Int32 kend) const
{
    if (MeiksinInitFailed()){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::GetIgmEndIndex: igm initialization not done");
    }

    Int32 Igm_kend = -1;
    // get last index in spectral axis where igm can be applied
    TAxisSampleList::const_iterator istart = m_SpectralAxis.GetSamplesVector().begin()+kstart;
    TAxisSampleList::const_iterator iend = m_SpectralAxis.GetSamplesVector().begin()+kend+1;
    TAxisSampleList::const_iterator it = std::upper_bound(istart, iend, m_igmCorrectionMeiksin->GetLambdaMax());
    return (it == istart) ? -1 : it - 1 - m_SpectralAxis.GetSamplesVector().begin(); // should be -1 if not applicable (lambdamax< lmabd[kstart]) 
}

void CTemplate::ScaleFluxAxis(Float64 amplitude){

    CSpectrum::ScaleFluxAxis(amplitude);
    if (CheckIsmIgmEnabled()) 
        m_NoIsmIgmFluxAxis *= amplitude;
}
void  CTemplate::GetIsmIgmIdxList(Int32 opt_extinction,
                            Int32 opt_dustFitting,
                            TInt32List& MeiksinList, //return 
                            TInt32List& EbmvList, //return
                            bool keepigmism,
                            Float64 FitEbmvCoeff,
                            Int32 FitMeiksinIdx)const
{
    if (MeiksinInitFailed() && opt_extinction){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::GetIsmIgmIdxList: missing Meiksin initialization");
    }
    if(CalzettiInitFailed() && opt_dustFitting != -1){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::GetIsmIgmIdxList: missing Calzetti initialization");
    }

    Int32 MeiksinListSize = 1;
    if(opt_extinction && !keepigmism)
    {
        MeiksinListSize = m_igmCorrectionMeiksin->GetIdxCount();
    }
    MeiksinList.resize(MeiksinListSize);
    if(opt_extinction)
    {
        if(keepigmism)
            MeiksinList[0] = FitMeiksinIdx;
        else
            std::iota(MeiksinList.begin(), MeiksinList.end(), 0);
    }else{
        MeiksinList[0] = -1;
    }
    Int32 EbmvListSize = 1;
    if(!keepigmism && (opt_dustFitting==-10 || opt_dustFitting!=-1))//==-10 for TF and !=-1 for PCA
    {
        EbmvListSize = m_ismCorrectionCalzetti->GetNPrecomputedEbmvCoeffs();
    }

    EbmvList.resize(EbmvListSize);
    if(opt_dustFitting!=-1){
        if(keepigmism)
            EbmvList[0] = m_ismCorrectionCalzetti->GetEbmvIndex(FitEbmvCoeff);
        else if(opt_dustFitting==-10 || !keepigmism)
                std::iota(EbmvList.begin(), EbmvList.end(), 0);
            else 
                EbmvList[0] = opt_dustFitting;
    }else{
        EbmvList[0] = -1;
    } 
}

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
#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_
#define _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h"
#include "RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h"
#include "RedshiftLibrary/log/log.h"
#include <stdexcept>
#include <string>
#include <map>
#include <iostream>
namespace NSEpic
{
class CTemplateCatalog;
class CTemplate : public CSpectrum
{

public:

    CTemplate() = default;
    CTemplate( const std::string& name, const std::string& category );
    CTemplate( const std::string& name, const std::string& category,
	       CSpectrumSpectralAxis spectralAxis, CSpectrumFluxAxis fluxAxis);
    CTemplate(const CTemplate& other);
    CTemplate(CTemplate&& other);
    CTemplate(const CTemplate& other, const TFloat64List & mask);
    CTemplate& operator=(const CTemplate& other); 
    CTemplate& operator=(CTemplate&& other); 
    ~CTemplate()=default;

    // override Flux Setters to reset ism/igm
    void SetFluxAxis(const CSpectrumFluxAxis & fluxaxis) override;
    void SetFluxAxis(CSpectrumFluxAxis && fluxaxis) override;
    void SetSpectralAndFluxAxes(CSpectrumSpectralAxis spcaxis, CSpectrumFluxAxis fluxaxis) override;

    // override spectral axis Setters to reset ism/igm (since depends on wavelength)
    void SetSpectralAxis(const CSpectrumSpectralAxis & spectralaxis) override;
    void SetSpectralAxis(CSpectrumSpectralAxis && spectralaxis) override;

    // override changing component type to reset ism/igm
    void SetType(const CSpectrum::EType type) override;
    void SetType(const CSpectrum::EType type) const override;

    const std::string&  GetCategory() const;

    Bool Save(const char *filePath ) const;

    bool ApplyDustCoeff(Int32 kDust);
    bool ApplyMeiksinCoeff(Int32 meiksinIdx); 
    void ScaleFluxAxis(Float64 amplitude);
    Int32 GetIsmCoeff() const;
    Int32 GetIgmCoeff() const;
    const TFloat64List& GetcomputedDustCoeffs() const;
    const TFloat64List& GetcomputedMeiksinCoeffs() const;

    void GetIsmIgmRangeIndex(Int32& begin, Int32& end) const;
    Int32 GetIgmEndIndex() const;

    void InitIsmIgmConfig( Float64 redshift,
                           const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti = nullptr,
                           const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin=nullptr);
    void InitIsmIgmConfig(const TFloat64Range& lbdaRange, Float64 redshift,
                           const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti = nullptr,
                           const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin=nullptr);
    void InitIsmIgmConfig(Int32 kstart, Int32 kend, Float64 redshift,
                           const std::shared_ptr<CSpectrumFluxCorrectionCalzetti>& ismCorrectionCalzetti = nullptr,
                           const std::shared_ptr<CSpectrumFluxCorrectionMeiksin>& igmCorrectionMeiksin=nullptr);
    void DisableIsmIgm();

    bool CheckIsmIgmEnabled() const {return !m_NoIsmIgmFluxAxis.isEmpty();};
    bool CalzettiInitFailed() const;
    bool MeiksinInitFailed() const;

    std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;
    std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;
    void   GetIsmIgmIdxList( Int32 opt_extinction,
                      Int32 opt_dustFitting,
                      TInt32List& MeiksinList, //return 
                      TInt32List& EbmvList, //return
                      Bool keepigmism = 0,
                      Float64 FitEbmvCoeff = NAN,
                      Int32 FitMeiksinIdx = -1) const;
private:

    std::string     m_Category;

    Int32   m_kDust = -1; 
    Int32   m_meiksinIdx = -1;
    Int32   m_meiksinRedshiftIdx = -1;

    Int32 m_IsmIgm_kstart = -1, m_Ism_kend = -1, m_Igm_kend = -1;
    CSpectrumFluxAxis   m_NoIsmIgmFluxAxis;

    //below vectors should be updated each time we change m_kDust, m_meiksinIdx for a specific redshift
    TFloat64List m_computedDustCoeff; //vector of spectrum size containing computed dust coeff at m_kDust and this for all lambdas in the spectrum
    TFloat64List m_computedMeiksingCoeff; //vector of spectrum size containing computed igm coeff at a specific Z at m_meiksin and this for all lambdas in the spectrum
};

// override Flux Setters to reset ism/igm
inline
void CTemplate::SetFluxAxis(const CSpectrumFluxAxis & fluxaxis)
{
    m_NoIsmIgmFluxAxis.clear();
    CSpectrum::SetFluxAxis(fluxaxis);
}

inline
void CTemplate::SetFluxAxis(CSpectrumFluxAxis && fluxaxis)
{
    m_NoIsmIgmFluxAxis.clear();
    CSpectrum::SetFluxAxis(std::move(fluxaxis));

}

inline
void CTemplate::SetSpectralAndFluxAxes(CSpectrumSpectralAxis spcaxis, CSpectrumFluxAxis fluxaxis)
{
    m_NoIsmIgmFluxAxis.clear(); 
    CSpectrum::SetSpectralAndFluxAxes(std::move(spcaxis), std::move(fluxaxis));
}

// override spectral axis Setters to reset ism/igm (since depends on wavelength)
inline
void CTemplate::SetSpectralAxis(const CSpectrumSpectralAxis & spectralaxis)
{
    m_NoIsmIgmFluxAxis.clear();
    CSpectrum::SetSpectralAxis(spectralaxis);
}

inline    
void CTemplate::SetSpectralAxis(CSpectrumSpectralAxis && spectralaxis)
{
    m_NoIsmIgmFluxAxis.clear();
    CSpectrum::SetSpectralAxis(std::move(spectralaxis));
}

// override changing component type to reset ism/igm
inline
void CTemplate::SetType(const CSpectrum::EType type)
{
    if(m_spcType != type)
    {   
        DisableIsmIgm();
        CSpectrum::SetType(type);
    }
}

inline
void CTemplate::SetType(const CSpectrum::EType type) const 
{
    if(m_spcType != type)
    {   
        if (!CheckIsmIgmEnabled())
            CSpectrum::SetType(type);
        else
        {
            throw GlobalException(INTERNAL_ERROR,"CTemplate::SetType: cannot change component type when ism/igm enabled on a const CTemplate");
        }   
    }
}



inline
void CTemplate::DisableIsmIgm() 
{
    // check if not already disabled
    if (!m_NoIsmIgmFluxAxis.isEmpty()){
        GetFluxAxis_() = m_NoIsmIgmFluxAxis;
        m_NoIsmIgmFluxAxis.clear();
    }
}

inline
Int32 CTemplate::GetIsmCoeff() const
{
    if (!CheckIsmIgmEnabled()){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::GetIsmCoeff:  ismigm initialization not done");
    }
    return m_kDust;
}

inline
Int32 CTemplate::GetIgmCoeff() const
{
    if (!CheckIsmIgmEnabled()){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::GetIgmCoeff:  ismigm initialization not done");
    }
    return m_meiksinIdx;
}

inline
void CTemplate::GetIsmIgmRangeIndex(Int32& begin, Int32& ismend) const
{
    if (!CheckIsmIgmEnabled()){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::GetIsmIgmRangeIndex:  ism initialization not done");
    }    
    begin = m_IsmIgm_kstart;
    ismend = m_Ism_kend;
}

inline
Int32 CTemplate::GetIgmEndIndex() const
{
    if (!CheckIsmIgmEnabled() || MeiksinInitFailed()){
        throw GlobalException(INTERNAL_ERROR,"CTemplate::GetIgmEndIndex: igm initialization not done");
    }
    return m_Igm_kend;
}

inline const TFloat64List& CTemplate::GetcomputedDustCoeffs() const{ return m_computedDustCoeff;}
inline const TFloat64List& CTemplate::GetcomputedMeiksinCoeffs() const{ return m_computedMeiksingCoeff;};

typedef std::vector< std::shared_ptr<CTemplate> >          TTemplateRefList;
typedef std::vector< std::shared_ptr< const CTemplate> >   TTemplateConstRefList;

typedef std::map< std::string, TTemplateRefList >          TTemplatesRefDict;
typedef std::map< std::string, TTemplateConstRefList >     TTemplatesConstRefDict;
}

#endif

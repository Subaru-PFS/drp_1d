#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_
#define _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h>
#include <RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h>
#include <RedshiftLibrary/log/log.h>
#include <stdexcept>
#include <string>
#include <map>

namespace NSEpic
{
class CTemplateCatalog;
class CTemplate : public CSpectrum
{

public:

    CTemplate();
    CTemplate( const std::string& name, const std::string& category );
    CTemplate( const std::string& name, const std::string& category,
	       CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis);
    CTemplate(const CTemplate& other);
    CTemplate& operator=(const CTemplate& other);
    ~CTemplate();

    const std::string&  GetCategory() const;
    const std::string&  GetName() const;
    Int32 GetTemplateByName(const CTemplateCatalog& tplCatalog,
                            const TStringList& tplCategoryList,
                            const std::string tplName,
                            CTemplate& tpl);
    Bool Save(const char *filePath ) const;

    const CSpectrumFluxAxis&        GetFluxAxisIsmIgm() const;
    CSpectrumFluxAxis&              GetFluxAxisIsmIgm();
    //Calzetti extinction
    bool ApplyDustCoeff(Int32 kDust);
    bool ApplyMeiksinCoeff(Int32 meiksinIdx, Float64 redshift); 
    Int32 GetIsmCoeff() const;
    Int32 GetIgmCoeff() const;
    void SetRequiredCorrections(Int32 nbcorrections);
    void DecrementCorrectionState();
    //void SetIsmIgmLambdaRange(Int32 kstart, Int32 kend) const;
    void SetFluxCorrectionIsmIgm(const std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti, 
                                 const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin);
    bool ReinitIsmIgmFlux();
private:

    std::string     m_Category;
    Int32   m_kDust = -1; //définie comme mutable pour pouvoir la changer dans Apply..coeff(), sinon ca ne marche pas
    Int32   m_meiksinIdx = -1;
    std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;
    //IGM meiksin
    std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;
    CSpectrumFluxAxis   m_FluxAxisIsmIgm;//flux on which is applied the igm and ism correction
    //below vectors should be updated each time we change m_kDust, m_meiksinIdx for a specific redshift
    TFloat64List m_computedDustCoeff; //vector of spectrum size containing computed dust coeff at m_kDust and this for all lambdas in the spectrum
    TFloat64List m_computedMeiksingCoeff; //vector of spectrum size containing computed igm coeff at a specific Z at m_meiksin and this for all lambdas in the spectrum
    
protected:
    Int32  m_IsmIgmApplied = -1;
    Int32  m_IsmIgmAppliedStatus = -1;
};
inline
const CSpectrumFluxAxis& CTemplate::GetFluxAxisIsmIgm() const
{
    if(m_kDust ==-1 && m_meiksinIdx == -1 && m_IsmIgmApplied == -1)
        return m_FluxAxis;
    else
        return m_FluxAxisIsmIgm;
}
inline
CSpectrumFluxAxis& CTemplate::GetFluxAxisIsmIgm()
{
    if(m_kDust ==-1 && m_meiksinIdx == -1 && m_IsmIgmApplied == -1)
        return m_FluxAxis;
    else
        return m_FluxAxisIsmIgm;
}
inline
Int32 CTemplate::GetIsmCoeff() const
{
    return m_kDust;
}
inline
Int32 CTemplate::GetIgmCoeff() const
{
    return m_meiksinIdx;
}

/**
 * @m_IsmIgmApplied
 * -1: no corrections applied; 
 * 0: all required corrections are applied;
 * 1: a first correction should be applied; 
 * 2: a second correction should be applied;
 * //each time a correction is applied we decrement m_IsmIgmApplied; 
 * once we reach 0== all corrections are applied, then we can multiply by the fluxaxis
*/
inline
void CTemplate::SetRequiredCorrections(Int32 nbcorrections ) 
{   
    m_IsmIgmApplied = nbcorrections; //max value
    m_IsmIgmAppliedStatus = nbcorrections;
    return;
}
inline
void CTemplate::DecrementCorrectionState()
{   
    m_IsmIgmAppliedStatus = m_IsmIgmAppliedStatus -1; 
    if(m_IsmIgmAppliedStatus < -1){
        Log.LogError("No corrections were setup in advance. CTemplate::SetRequiredCorrections should be called first! Aborting");
        //throw runtime_error("No ISM or IGM corrections were setup. CTemplate::SetRequiredCorrections should be called first! Aborting");
    }
    return;
}
typedef std::vector< std::shared_ptr<CTemplate> >          TTemplateRefList;
typedef std::vector< std::shared_ptr< const CTemplate> >     TTemplateConstRefList;

typedef std::map< std::string, TTemplateRefList >          TTemplatesRefDict;
typedef std::map< std::string, TTemplateConstRefList >     TTemplatesConstRefDict;
}

#endif

#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_
#define _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/fluxcorrectionmeiksin.h>
#include <RedshiftLibrary/spectrum/fluxcorrectioncalzetti.h>
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
    bool ApplyDustCoeff(Int32 kDust) const;
    bool ApplyMeiksinCoeff(Int32 meiksinIdx, Float64 redshift) const; 
    Int32 GetIsmCoeff() const;
    Int32 GetIgmCoeff() const;
    //void SetIsmIgmLambdaRange(Int32 kstart, Int32 kend) const;
    void SetFluxCorrectionIsmIgm(const std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti, 
                                 const std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin) const;
    bool ReinitIsmIgmFlux();
private:

    std::string     m_Category;
    mutable Int32   m_kDust = -1; //définie comme mutable pour pouvoir la changer dans Apply..coeff(), sinon ca ne marche pas
    mutable Int32   m_meiksinIdx = -1;
    mutable std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;
    //IGM meiksin
    mutable std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;
    mutable CSpectrumFluxAxis   m_FluxAxisIsmIgm;//flux on which is applied the igm and ism correction
    //below vectors should be updated each time we change m_kDust, m_meiksinIdx for a specific redshift
    mutable TFloat64List computedDustCoeff; //vector of spectrum size containing computed dust coeff at m_kDust and this for all lambdas in the spectrum
    mutable TFloat64List computedMeiksingCoeff; //vector of spectrum size containing computed igm coeff at a specific Z at m_meiksin and this for all lambdas in the spectrum
    
protected:
    mutable Bool   m_IsmIgmApplied = false;//not yet used
};
inline
const CSpectrumFluxAxis& CTemplate::GetFluxAxisIsmIgm() const
{
    if(m_kDust ==-1 && m_meiksinIdx == -1)
        return m_FluxAxis;
    else
        return m_FluxAxisIsmIgm;
}
inline
CSpectrumFluxAxis& CTemplate::GetFluxAxisIsmIgm()
{
    if(m_kDust ==-1 && m_meiksinIdx == -1)
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
typedef std::vector< std::shared_ptr<CTemplate> >          TTemplateRefList;
typedef std::vector< std::shared_ptr< const CTemplate> >     TTemplateConstRefList;

typedef std::map< std::string, TTemplateRefList >          TTemplatesRefDict;
typedef std::map< std::string, TTemplateConstRefList >     TTemplatesConstRefDict;
}

#endif

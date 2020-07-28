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
    void SetIsmIgmLambdaRange(Int32 kstart, Int32 kend) const;
    void SetFluxCorrectionIsmIgm(std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti, std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin) const;
private:

    std::string     m_Category;
    std::string     m_Name;
    mutable Int32   m_kstart;
    mutable Int32   m_kend;
    mutable Int32   m_kDust; //définie comme mutable pour pourvoir la changer dans Apply..coeff(), sinon ca ne marche pas
    mutable Int32   m_meiksinIdx;
    mutable std::shared_ptr<CSpectrumFluxCorrectionCalzetti> m_ismCorrectionCalzetti;
    //IGM meiksin
    mutable std::shared_ptr<CSpectrumFluxCorrectionMeiksin> m_igmCorrectionMeiksin;
    CSpectrumFluxAxis   m_FluxAxisIsmIgm;//flux on which is applied the igm and ism correction
};
inline
const CSpectrumFluxAxis& CTemplate::GetFluxAxisIsmIgm() const
{
    return m_FluxAxisIsmIgm;
}
inline
CSpectrumFluxAxis& CTemplate::GetFluxAxisIsmIgm()
{
    return m_FluxAxisIsmIgm;
}

typedef std::vector< std::shared_ptr<CTemplate> >          TTemplateRefList;
typedef std::vector< std::shared_ptr< const CTemplate> >     TTemplateConstRefList;

typedef std::map< std::string, TTemplateRefList >          TTemplatesRefDict;
typedef std::map< std::string, TTemplateConstRefList >     TTemplatesConstRefDict;
}

#endif

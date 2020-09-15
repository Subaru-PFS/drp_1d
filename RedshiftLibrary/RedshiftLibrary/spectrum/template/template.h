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
#include <iostream>
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
    const CSpectrumFluxAxis&    GetFluxAxis() const;
    CSpectrumFluxAxis&          GetFluxAxis();
    const CSpectrumFluxAxis&    GetFluxAxisWithoutIsmIgm() const;
    CSpectrumFluxAxis&          GetFluxAxisWithoutIsmIgm();

    Bool Save(const char *filePath ) const;

    //Calzetti extinction
    bool ApplyDustCoeff(Int32 kDust);
    bool ApplyMeiksinCoeff(Int32 meiksinIdx, Float64 redshift); 

    Int32 GetIsmCoeff() const;
    Int32 GetIgmCoeff() const;

    void SetIsmIgmLambdaRange(Int32 kstart, Int32 kend);
    bool InitIsmIgmConfig();

    CSpectrumFluxCorrectionCalzetti m_ismCorrectionCalzetti;
    CSpectrumFluxCorrectionMeiksin m_igmCorrectionMeiksin;
private:

    std::string     m_Category;
    Int32   m_kDust = -1; //définie comme mutable pour pouvoir la changer dans Apply..coeff(), sinon ca ne marche pas
    Int32   m_meiksinIdx = -1;
    Float64 m_redshiftMeiksin = -1;

    Int32 m_kstart = -1, m_kend = -1;
    CSpectrumFluxAxis   m_FluxAxisIsmIgm;//flux on which is applied the igm and ism correction
    //below vectors should be updated each time we change m_kDust, m_meiksinIdx for a specific redshift
    TFloat64List m_computedDustCoeff; //vector of spectrum size containing computed dust coeff at m_kDust and this for all lambdas in the spectrum
    TFloat64List m_computedMeiksingCoeff; //vector of spectrum size containing computed igm coeff at a specific Z at m_meiksin and this for all lambdas in the spectrum
    
protected:
    Int32  m_IsmIgmApplied = -1;
};

//override spectrum flux getters to return the corrected flux rather than the raw flux
//this is required when saving the spectrum model corresponding to a template. By default
//it is considered as a spectrum..thus we lose the m_FluxAxisIsmIgm variable.
inline
const CSpectrumFluxAxis& CTemplate::GetFluxAxis() const
{
    if(m_kDust ==-1 && m_meiksinIdx == -1)
        return m_FluxAxis;
    else{
        return m_FluxAxisIsmIgm;
    }
}

inline
CSpectrumFluxAxis& CTemplate::GetFluxAxis()
{
    if(m_kDust ==-1 && m_meiksinIdx == -1){
        return m_FluxAxis;
    }
    else{
        return m_FluxAxisIsmIgm;
    } 
}

inline
const CSpectrumFluxAxis& CTemplate::GetFluxAxisWithoutIsmIgm() const
{
    return m_FluxAxis;
   
}
inline
CSpectrumFluxAxis& CTemplate::GetFluxAxisWithoutIsmIgm()
{
    return m_FluxAxis;
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

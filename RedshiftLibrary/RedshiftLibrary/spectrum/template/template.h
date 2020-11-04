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

    const CSpectrumFluxAxis&    GetFluxAxisWithoutIsmIgm() const;
    CSpectrumFluxAxis&          GetFluxAxisWithoutIsmIgm();

    Bool Save(const char *filePath ) const;

    bool ApplyDustCoeff(Int32 kDust);
    bool ApplyMeiksinCoeff(Int32 meiksinIdx, Float64 redshift); 
    void ScaleFluxAxis(Float64 amplitude);
    Int32 GetIsmCoeff() const;
    Int32 GetIgmCoeff() const;

    void SetIsmIgmLambdaRange(TFloat64Range& lbdaRange);
    void GetIsmIgmRangeIndex(Int32& begin, Int32& end);
    bool InitIsmIgmConfig();

    CSpectrumFluxCorrectionCalzetti m_ismCorrectionCalzetti;
    CSpectrumFluxCorrectionMeiksin m_igmCorrectionMeiksin;
private:

    std::string     m_Category;

    Int32   m_kDust = -1; 
    Int32   m_meiksinIdx = -1;
    Float64 m_redshiftMeiksin = -1;

    Int32 m_IsmIgm_kstart = -1, m_IsmIgm_kend = -1;
    CSpectrumFluxAxis   m_NoIsmIgmFluxAxis;
    //below vectors should be updated each time we change m_kDust, m_meiksinIdx for a specific redshift
    TFloat64List m_computedDustCoeff; //vector of spectrum size containing computed dust coeff at m_kDust and this for all lambdas in the spectrum
    TFloat64List m_computedMeiksingCoeff; //vector of spectrum size containing computed igm coeff at a specific Z at m_meiksin and this for all lambdas in the spectrum
};


inline
const CSpectrumFluxAxis& CTemplate::GetFluxAxisWithoutIsmIgm() const
{
    return m_NoIsmIgmFluxAxis;
   
}
inline
CSpectrumFluxAxis& CTemplate::GetFluxAxisWithoutIsmIgm()
{
    return m_NoIsmIgmFluxAxis;
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
inline
void CTemplate::GetIsmIgmRangeIndex(Int32& begin, Int32& end){
    begin = m_IsmIgm_kstart;
    end   = m_IsmIgm_kend;
}
typedef std::vector< std::shared_ptr<CTemplate> >          TTemplateRefList;
typedef std::vector< std::shared_ptr< const CTemplate> >   TTemplateConstRefList;

typedef std::map< std::string, TTemplateRefList >          TTemplatesRefDict;
typedef std::map< std::string, TTemplateConstRefList >     TTemplatesConstRefDict;
}

#endif

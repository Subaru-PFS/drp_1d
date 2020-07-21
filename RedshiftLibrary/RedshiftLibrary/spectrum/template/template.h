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
    Bool Save(const char *filePath ) const;
    Int32  RebinTemplate( const CSpectrum& spectrum, 
                                Float64 redshift,
                                const TFloat64Range& lambdaRange,
                                std::string opt_interp,
                                CSpectrumSpectralAxis& spcSpectralAxis_restframe,
                                CSpectrum& itplTplSpectrum,
                                CMask& itplMask,
                                TFloat64Range& currentRange,
                                Float64& overlaprate,
                                Float64 overlapThreshold,
                                TFloat64List& YtplRawBuffer) const;

    //Calzetti extinction
    bool ApplyDustCoeff(Int32 kDust, const TAxisSampleList & Xtpl, TAxisSampleList & Ytpl/*, TFloat64List Ytpl*/, Int32 kstart, Int32 kend, CSpectrumFluxCorrectionCalzetti* ismCorrectionCalzetti) const;
    
    bool ApplyMeiksinCoeff(Int32 meiksinIdx, const TAxisSampleList & Xtpl, TAxisSampleList & Ytpl/*, TFloat64List Ytpl*/, Int32 kstart, Int32 kend, Int32 redshiftIdx, CSpectrumFluxCorrectionMeiksin* igmCorrectionMeiksin) const;

   
private:

    std::string     m_Category;
    std::string     m_Name;
};

typedef std::vector< std::shared_ptr<CTemplate> >          TTemplateRefList;
typedef std::vector< std::shared_ptr< const CTemplate> >     TTemplateConstRefList;

typedef std::map< std::string, TTemplateRefList >          TTemplatesRefDict;
typedef std::map< std::string, TTemplateConstRefList >     TTemplatesConstRefDict;
}

#endif

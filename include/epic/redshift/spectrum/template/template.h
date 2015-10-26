#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_
#define _REDSHIFT_SPECTRUM_TEMPLATE_TEMPLATE_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/ref.h>
#include <epic/core/common/constref.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <string>
#include <map>

namespace NSEpic
{

class CTemplate : public CSpectrum
{

public:

    CTemplate();
    CTemplate( const std::string& name, const std::string& category );
    ~CTemplate();

    const std::string&  GetCategory() const;
    const std::string&  GetName() const;

private:

    std::string     m_Category;
    std::string     m_Name;
};

typedef std::vector< CRef<CTemplate> >          TTemplateRefList;
typedef std::vector< CConstRef<CTemplate> >     TTemplateConstRefList;

typedef std::map< std::string, TTemplateRefList >          TTemplatesRefDict;
typedef std::map< std::string, TTemplateConstRefList >     TTemplatesConstRefDict;
}

#endif

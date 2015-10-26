#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_
#define _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>
#include <epic/core/common/ref.h>
#include <epic/redshift/spectrum/template/template.h>

#include <boost/filesystem.hpp>
#include <vector>

namespace NSEpic
{

class CTemplateCatalog : public CManagedObject
{

    DEFINE_MANAGED_OBJECT( CTemplateCatalog )

public:

    CTemplateCatalog();
    ~CTemplateCatalog();

    Bool Add( CTemplate& r );
    Bool Add( const char* templatePath, const std::string& category );
    Bool Load( const char* filePath );

    const CTemplate& GetTemplate( const std::string& category, UInt32 i ) const;
    const CTemplate& GetTemplateWithoutContinuum( const std::string& category, UInt32 i ) const;

    TTemplateRefList GetTemplate( const TStringList& categoryList ) const;
    TTemplateRefList GetTemplateWithoutContinuum(  const TStringList& categoryList  ) const;

    TStringList GetCategoryList() const;

    UInt32 GetTemplateCount( const std::string& category ) const;

private:

    Bool                    LoadCategory( const boost::filesystem::path& dirPath, const std::string& category );

    TTemplatesRefDict        m_List;
    TTemplatesRefDict        m_ListWithoutCont;
};


}

#endif

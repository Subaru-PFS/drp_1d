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
    Bool Add( const char* templatePath, CTemplate::ECategory category );
    Bool Load( const char* filePath );

    const CTemplate& GetTemplate( CTemplate::ECategory category, UInt32 i ) const;
    const CTemplate& GetTemplateWithoutContinuum( CTemplate::ECategory category, UInt32 i ) const;

    TTemplateRefList GetTemplate( const TTemplateCategoryList& categoryList ) const;
    TTemplateRefList GetTemplateWithoutContinuum(  const TTemplateCategoryList& categoryList  ) const;

    TTemplateCategoryList GetCategoryList() const;

    UInt32 GetTemplateCount( CTemplate::ECategory category ) const;

private:

    Bool                    LoadCategory( const boost::filesystem::path& dirPath, CTemplate::ECategory category );
    CTemplate::ECategory    ConvertStringToCategory( const std::string& category );
    TTemplateRefList        m_List[CTemplate::nCategory_Count];
    TTemplateRefList        m_ListWithoutCont[CTemplate::nCategory_Count];
};


}

#endif

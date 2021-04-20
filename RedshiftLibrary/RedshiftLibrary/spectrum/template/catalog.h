#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_
#define _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/template/template.h>

#include <boost/filesystem.hpp>
#include <vector>

namespace NSEpic
{

class CTemplateCatalog
{

public:

    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8.0, std::string waveletsDFBinPath="" );
    ~CTemplateCatalog();

    void Add( std::shared_ptr<CTemplate> );
    void Add( const char* templatePath, const std::string& category );
    void Load( const char* filePath );
    Bool Save(const char* filePath , Bool saveWithoutContinuum=true);

    const CTemplate& GetTemplate( const std::string& category, UInt32 i ) const;
    const CTemplate& GetTemplateByName(const TStringList& tplCategoryList, const std::string tplName ) const;

    TTemplateConstRefList GetTemplate( const TStringList& categoryList ) const;
    TTemplateRefList GetTemplate( const TStringList& categoryList );

    TStringList GetCategoryList() const;

    UInt32 GetTemplateCount( const std::string& category ) const;
    
    static TTemplateConstRefList const_TTemplateRefList_cast(const TTemplateRefList & list);

private:
    // this const version must stay private, since it returns non const templates.
    TTemplateRefList GetTemplate_( const TStringList& categoryList ) const; 

    Bool                     LoadCategory( const boost::filesystem::path& dirPath, const std::string& category );

    TTemplatesRefDict        m_List;

    std::string m_continuumRemovalMethod;
    Float64 m_continuumRemovalMedianKernelWidth;
    Float64 m_continuumRemovalWaveletsNScales;
    std::string m_continuumRemovalWaveletsBinPath;

};

/**
 * Returns the contents of the i-th entry in the category item of m_List.
 */
inline const CTemplate& CTemplateCatalog::GetTemplate( const std::string& category, UInt32 i ) const
{
    return *m_List.at( category )[i];
}

// non const getter returning mutable templates
inline 
TTemplateRefList CTemplateCatalog::GetTemplate( const TStringList& categoryList )
{
    return GetTemplate_(categoryList);
}

//  const getter returning const templates
inline
TTemplateConstRefList CTemplateCatalog::GetTemplate( const TStringList& categoryList ) const
{
    return const_TTemplateRefList_cast( GetTemplate_(categoryList)); 
}

}

#endif

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

    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8.0, std::string waveletsDFBinPath="", Bool scale = 0 );
    
    void                    SetCurrentScale(std::string scale);
    void                    Add( std::shared_ptr<CTemplate> , std::string scale ="lin");
    const CTemplate&        GetTemplate( const std::string& category, UInt32 i ) const;
    const CTemplate&        GetTemplateByName(const TStringList& tplCategoryList, const std::string tplName ) const;

    TTemplateConstRefList GetTemplate( const TStringList& categoryList ) const;
    TTemplateRefList GetTemplate( const TStringList& categoryList );

    TStringList GetCategoryList() const;

    UInt32 GetTemplateCount( const std::string& category ) const;
    
    static TTemplateConstRefList const_TTemplateRefList_cast(const TTemplateRefList & list);

    TStringList             GetCategoryList() const;
    UInt32                  GetTemplateCount( const std::string& category ) const;
    void                    InitIsmIgm(const std::string & calibrationPath);
private:
    // this const version must stay private, since it returns non const templates.
    TTemplateRefList GetTemplate_( const TStringList& categoryList ) const; 

    //Bool                    LoadCategory( const boost::filesystem::path& dirPath, const std::string& category );

    TTemplatesRefDict        m_List;
    TTemplatesRefDict        m_ListRebinned;
    //we need to find a mechanism to select  m_list vs m_listRebinned 
    Bool                     m_logscale = 0;//non-log by default

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
/* //below functions aim at avoid using if..else to access the right categoryList
inline 
TTemplateRefList  CTemplateCatalog::GetCategoryList(){
    if(!m_logscale)
        return *m_List.at( category )[i];
    else
        return m_ListRebinned;       
}
inline 
TTemplateRefList  CTemplateCatalog::GetCategoryList(Bool logscale){
    if(!logscale)
        return m_List;
    else
        return m_ListRebinned;       
}
}

#endif

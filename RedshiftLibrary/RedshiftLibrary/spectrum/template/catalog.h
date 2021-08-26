#ifndef _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_
#define _REDSHIFT_SPECTRUM_TEMPLATE_CATALOG_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/template/template.h"

#include <boost/filesystem.hpp>
#include <vector>

namespace NSEpic
{

class CTemplateCatalog
{

public:

    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8.0, std::string waveletsDFBinPath="", Bool sampling = 0 );
    
    void                    Add( const std::shared_ptr<CTemplate> & tpl);

    std::shared_ptr<const CTemplate>        GetTemplate( const std::string& category, UInt32 i ) const;
    std::shared_ptr<const CTemplate>        GetTemplateByName(const TStringList& tplCategoryList, const std::string tplName ) const;
    void                                    SetTemplate( const std::shared_ptr<CTemplate> & tpl,  UInt32 i);

    TTemplateConstRefList GetTemplateList( const TStringList& categoryList ) const;
    TTemplateRefList GetTemplateList( const TStringList& categoryList );
    
    static TTemplateConstRefList const_TTemplateRefList_cast(const TTemplateRefList & list);

    TStringList             GetCategoryList() const;
    UInt32                  GetTemplateCount( const std::string& category ) const;
    void                    InitIsmIgm(const std::string & calibrationPath, 
                                       std::shared_ptr<const CParameterStore> parameterStore,
                                       const std::shared_ptr<const CLSF>& lsf);
    mutable Bool            m_logsampling = 0;//non-log by default

private:
    // this const version must stay private, since it returns non const templates.
    TTemplateRefList GetTemplateList_( const TStringList& categoryList ) const; 

    //Bool                    LoadCategory( const boost::filesystem::path& dirPath, const std::string& category );
          TTemplatesRefDict &    GetList();//using m_sampling
    const TTemplatesRefDict &    GetList() const;

    TTemplatesRefDict        m_List;
    TTemplatesRefDict        m_ListRebinned;

    std::string m_continuumRemovalMethod;
    Float64 m_continuumRemovalMedianKernelWidth;
    Float64 m_continuumRemovalWaveletsNScales;
    std::string m_continuumRemovalWaveletsBinPath;

};

/**
 * Returns the contents of the i-th entry in the category item of m_List.
 */
inline 
std::shared_ptr<const CTemplate> CTemplateCatalog::GetTemplate( const std::string& category, UInt32 i ) const
{
    return GetList().at( category )[i];   
}

// non const getter returning mutable templates
inline 
TTemplateRefList CTemplateCatalog::GetTemplateList( const TStringList& categoryList )
{
    return GetTemplateList_(categoryList);
}

//  const getter returning const templates
inline
TTemplateConstRefList CTemplateCatalog::GetTemplateList( const TStringList& categoryList ) const
{
    return const_TTemplateRefList_cast( GetTemplateList_(categoryList)); 
}


 //below functions aim at avoid using if..else to access the right categoryList
inline 
const TTemplatesRefDict & CTemplateCatalog::GetList() const
{
    return const_cast<CTemplateCatalog*>(this)->GetList();
}

inline 
TTemplatesRefDict & CTemplateCatalog::GetList()
{
    return m_logsampling? m_ListRebinned : m_List;
}

}

#endif

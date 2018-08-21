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

    CTemplateCatalog( std::string cremovalmethod="Median", Float64 mediankernelsize=75.0, Float64 waveletsScales=8, std::string waveletsDFBinPath="");
    ~CTemplateCatalog();

    void Add( std::shared_ptr<CTemplate> );
    void Add( const char* templatePath, const std::string& category );
    void Load( const char* filePath );
    Bool Save(const char* filePath , Bool saveWithoutContinuum=true);

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


    std::string m_continuumRemovalMethod;
    Float64 m_continuumRemovalMedianKernelWidth;
    Float64 m_continuumRemovalWaveletsNScales;
    std::string m_continuumRemovalWaveletsBinPath;

};


}

#endif

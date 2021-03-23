#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/spectrum/io/genericreader.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/continuum/median.h>
#include <RedshiftLibrary/continuum/irregularsamplingmedian.h>
#include <RedshiftLibrary/continuum/waveletsdf.h>

#include <RedshiftLibrary/log/log.h>

#include <boost/filesystem.hpp>

#include <string>

using namespace NSEpic;
using namespace std;
using namespace boost::filesystem;


/**
 * Variable instantiator constructor.
 */
CTemplateCatalog::CTemplateCatalog( std::string cremovalmethod, Float64 mediankernelsize, Float64 waveletsScales, std::string waveletsDFBinPath, Bool scale )
{
    m_continuumRemovalMethod = cremovalmethod;
    m_continuumRemovalMedianKernelWidth = mediankernelsize;
    m_continuumRemovalWaveletsNScales = waveletsScales;
    m_continuumRemovalWaveletsBinPath = waveletsDFBinPath;
    m_logscale = scale;
}


TTemplateConstRefList CTemplateCatalog::const_TTemplateRefList_cast(const TTemplateRefList & list)
{
    TTemplateConstRefList const_list;
    for (auto tpl : list) const_list.push_back(tpl);

    return const_list;
}

/**
 * Returns a list containing all templates as enumerated in the categoryList input.
 */
TTemplateRefList CTemplateCatalog::GetTemplate_( const TStringList& categoryList ) const
{
    TTemplateRefList list;

    for( Int32 i=0; i<categoryList.size(); i++ )
    {
        for ( Int32 j=0; j<GetTemplateCount( categoryList[i] ); j++ )
        {
            if(!m_logscale)
                list.push_back( m_List.at( categoryList[i] )[j] );
            else
                list.push_back( m_ListRebinned.at( categoryList[i] )[j] );
        }
    }

    return list;
}

std::shared_ptr<const CTemplate>  CTemplateCatalog::GetTemplateByName(const TStringList& tplCategoryList, const std::string tplName ) const
{
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        for( UInt32 j=0; j<GetTemplateCount( tplCategoryList[i] ); j++ )
        {
            std::shared_ptr<const CTemplate> tpl = GetTemplate( tplCategoryList[i], j );
            if(tpl->GetName() == tplName){
                return tpl;
            }
        }
    }
    throw std::runtime_error("Could not find template with name");
}


/**
 * Get a list of strings with the contents of m_List.
 */
TStringList CTemplateCatalog::GetCategoryList() const
{
    TStringList l;
    
    if(!m_logscale)
        for( auto it = m_List.begin(); it != m_List.end(); it++ ) {
            l.push_back( it->first );
        }
    else
    {
        for( auto it = m_ListRebinned.begin(); it != m_List.end(); it++ ) {
            l.push_back( it->first );
        }
    }

    return l;
}

/**
 * Returns the size of the category entry in m_List.
 */
UInt32 CTemplateCatalog::GetTemplateCount( const std::string& category ) const
{   
    UInt32 l; 
    if(!m_logscale){
        if( m_List.find( category ) == m_List.end() )
            return 0;
        l = m_List.at( category ).size();
    }else{
        if( m_ListRebinned.find( category ) == m_List.end() )
            return 0;
        l = m_ListRebinned.at( category ).size();
    }
    return l;
}

/**
 * Adds the input to the list of templates, under its category. If the input doesn't have a category, function returns false. Also computes the template without continuum and adds it to the list of templates without continuum. Returns true.
 */
void CTemplateCatalog::Add( std::shared_ptr<CTemplate> r, std::string scale)
{
    if( r->GetCategory().empty() )
      throw runtime_error("Template has no category");
    if(scale == "lin")
        m_List[r->GetCategory()].push_back( r );
    if(scale == "log")
        m_ListRebinned[r->GetCategory()].push_back( r );
}

void CTemplateCatalog::InitIsmIgm(const std::string & calibrationPath)
{
    //ISM
    auto ismCorrectionCalzetti = std::make_shared<CSpectrumFluxCorrectionCalzetti>();
    ismCorrectionCalzetti->Init(calibrationPath, 0.0, 0.1, 10);
    //IGM
    auto igmCorrectionMeiksin = std::make_shared<CSpectrumFluxCorrectionMeiksin>();
    igmCorrectionMeiksin->Init(calibrationPath);

    //push in all templates
    for(std::string s : GetCategoryList()){ 
        TTemplateRefList  TplList = GetTemplate(TStringList{s});
        for (auto tpl : TplList)
        {
            tpl->m_ismCorrectionCalzetti = ismCorrectionCalzetti;
            if(s!="star")//no igm for stars
                tpl->m_igmCorrectionMeiksin = igmCorrectionMeiksin;
        }   
    }
}
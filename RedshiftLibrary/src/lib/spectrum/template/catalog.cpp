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
CTemplateCatalog::CTemplateCatalog( std::string cremovalmethod, Float64 mediankernelsize, Float64 waveletsScales, std::string waveletsDFBinPath )
{
    m_continuumRemovalMethod = cremovalmethod;
    m_continuumRemovalMedianKernelWidth = mediankernelsize;
    m_continuumRemovalWaveletsNScales = waveletsScales;
    m_continuumRemovalWaveletsBinPath = waveletsDFBinPath;
}

/**
 * Empty destructor.
 */
CTemplateCatalog::~CTemplateCatalog()
{

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
            list.push_back( m_List.at( categoryList[i] )[j] );
        }
    }

    return list;
}


const CTemplate&  CTemplateCatalog::GetTemplateByName(const TStringList& tplCategoryList, const std::string tplName ) const
{
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        for( UInt32 j=0; j<GetTemplateCount( tplCategoryList[i] ); j++ )
        {
            const CTemplate& tpl = GetTemplate( tplCategoryList[i], j );
            if(tpl.GetName() == tplName){
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

    for( auto it = m_List.begin(); it != m_List.end(); it++ ) {
        l.push_back( it->first );
    }

    return l;
}

/**
 * Returns the size of the category entry in m_List.
 */
UInt32 CTemplateCatalog::GetTemplateCount( const std::string& category ) const
{
    if( m_List.find( category ) == m_List.end() )
        return 0;

    return m_List.at( category ).size();
}

/**
 * Adds the input to the list of templates, under its category. If the input doesn't have a category, function returns false. Also computes the template without continuum and adds it to the list of templates without continuum. Returns true.
 */
void CTemplateCatalog::Add( std::shared_ptr<CTemplate> r )
{
    if( r->GetCategory().empty() )
      throw runtime_error("Template has no category");

    m_List[r->GetCategory()].push_back( r );
}

/**
 * Adds the templates in the input path, under the input category.
 */
//TODO [ml] templatePath should be a string or boost::path
void CTemplateCatalog::Add( const char* templatePath, const std::string& category )
{
    if ( !exists( templatePath ) )
      {
	Log.LogError("Path %s does not exist.", templatePath);
	throw std::runtime_error("Template path does not exist");
      }

    path p( templatePath );
    std::string name = p.leaf().generic_string();

    std::shared_ptr<CTemplate> tmpl = std::make_shared<CTemplate>(name, category);
    tmpl->SetContinuumEstimationMethod(m_continuumRemovalMethod);
    tmpl->SetMedianWinsize(m_continuumRemovalMedianKernelWidth);
    tmpl->SetDecompScales(m_continuumRemovalWaveletsNScales);
    tmpl->SetWaveletsDFBinPath(m_continuumRemovalWaveletsBinPath);

    CSpectrumIOGenericReader asciiReader;
    asciiReader.Read( templatePath, *tmpl );

    Add( tmpl );
}

/**
 * Loads input directory as a collection of categories of templates.
 */
void CTemplateCatalog::Load( const char* dirPath )
{
    path pathFound;

    if ( !exists( dirPath ) )
      {
	Log.LogError("Path %s does not exist.", dirPath);
	throw std::runtime_error("Template catalog path does not exist");
      }

    //TODO [ml] initialize end_itr (how can this work ???)
    directory_iterator end_itr;
    for ( directory_iterator itr( dirPath ); itr != end_itr; ++itr )
    {
        if ( is_directory( itr->status() ) )
        {
            path::iterator it = itr->path().end();
            it--;

            std::string category = (*it).generic_string();
            LoadCategory( itr->path(), category );
        }
    }
}

/**
 * Save catalog in the dirPath.
 */
Bool CTemplateCatalog::Save( const char* dirPath, Bool saveWithoutContinuum )
{
    if ( !exists( dirPath ) )
      {
	Log.LogError ( "Path %s does not exist.", dirPath );
        return false;
      }

    path dirBase (dirPath);
    TStringList tplCategoryList = GetCategoryList();
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];
        path dirCategory (category);
        path dirCategoryFull = dirBase / dirCategory;
        create_directories(dirCategoryFull);

        for( UInt32 j=0; j<GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = GetTemplate( category, j );

            std::string filePath = tpl.GetName();
            path file (filePath.c_str());
            path full_path = dirCategoryFull / file;
            if(saveWithoutContinuum)
            {
                CSpectrum::EType savetype = tpl.GetType();
                tpl.SetType(CSpectrum::nType_noContinuum);
                tpl.Save( full_path.c_str() );
                tpl.SetType(savetype);
            }
            else
            {
                tpl.Save( full_path.c_str() );
            }
        }
    }
    return true;
}

/**
 * Load files in dirPath under the input category.
 */
Bool CTemplateCatalog::LoadCategory( const path& dirPath, const std::string&  category )
{
    //TODO [ml] initialize end_itr (how can this work ???)
    directory_iterator end_itr;
    for ( directory_iterator itr( dirPath ); itr != end_itr; ++itr )
    {
        //hack limit the number of templates loaded for this category
        //Int32 ntplMax = 7;
        //if(GetTemplateCount(category)>=ntplMax)
        //{
        //    return true;
        //}

        //hack filter the template catalog by name
        //std::string name = itr->path().c_str();
        //std::size_t foundstra = name.find("Rebinned-NEW-E-extendeddataExtensionData");
        //if (foundstra==std::string::npos)
        //{
        //    continue;
        //}

        if ( !is_directory( itr->status() ) )
	  {
	    Add( itr->path().c_str(), category );
    }
    }
    Log.LogInfo ( "Loaded %d templates for category %s.", GetTemplateCount(category), category.c_str() );

    return true;
}

#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/spectrum/io/genericreader.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/continuum/median.h>
#include <epic/redshift/continuum/irregularsamplingmedian.h>

#include <epic/core/log/log.h>

#include <boost/filesystem.hpp>

#include <string>

using namespace NSEpic;
using namespace std;
using namespace boost::filesystem;


/**
 * Variable instantiator constructor.
 */
CTemplateCatalog::CTemplateCatalog( string cremovalmethod, Float64 mediankernelsize )
{
    m_continuumRemovalMethod = cremovalmethod;
    m_continuumRemovalMedianKernelWidth = mediankernelsize;
}

/**
 * Empty destructor.
 */
CTemplateCatalog::~CTemplateCatalog()
{

}

/**
 * Returns the contents of the i-th entry in the category item of m_List.
 */
const CTemplate& CTemplateCatalog::GetTemplate( const std::string& category, UInt32 i ) const
{
    return *m_List.at( category )[i];
}

/**
 * Returns the contents of the i-th entry in the category item of m_ListWithoutCont.
 */
const CTemplate& CTemplateCatalog::GetTemplateWithoutContinuum( const std::string& category, UInt32 i ) const
{
    return *m_ListWithoutCont.at( category )[i];
}

/**
 * Returns a list containing all templates as enumerated in the categoryList input.
 */
TTemplateRefList CTemplateCatalog::GetTemplate( const TStringList& categoryList ) const
{
    TTemplateRefList list;

    for( Int32 i=0; i<categoryList.size(); i++ )
    {
        for ( Int32 j=0; j<GetTemplateCount( categoryList[i] ); i++ )
        {
            list.push_back( m_List.at( categoryList[i] )[j] );
        }
    }

    return list;
}

/**
 * Returns a list containing all templates without continuum as enumerated in the categoryList input.
 */
TTemplateRefList CTemplateCatalog::GetTemplateWithoutContinuum( const TStringList& categoryList ) const
{
    TTemplateRefList list;

    for( Int32 i=0; i<categoryList.size(); i++ )
    {
        for ( Int32 j=0; j<GetTemplateCount( categoryList[i] ); i++ )
        {
            list.push_back( m_ListWithoutCont.at( categoryList[i] )[j] );
        }
    }

    return list;
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
Bool CTemplateCatalog::Add( std::shared_ptr<CTemplate> r )
{
    if( r->GetCategory().empty() )
        return false;

    m_List[r->GetCategory()].push_back( r );

    // Compute continuum substracted spectrum
    std::shared_ptr<CTemplate> tmplWithoutCont = std::shared_ptr<CTemplate>( new CTemplate( r->GetName().c_str(), r->GetCategory() ) );

    *tmplWithoutCont = *r;

    if( m_continuumRemovalMethod == "Median" )
      {
        CContinuumMedian continuum;
        continuum.SetMedianKernelWidth( m_continuumRemovalMedianKernelWidth );
        tmplWithoutCont->RemoveContinuum( continuum );
      }
    else
      {
        CContinuumIrregularSamplingMedian continuum;
        continuum.SetMedianKernelWidth( m_continuumRemovalMedianKernelWidth );
        tmplWithoutCont->RemoveContinuum( continuum );
      }

    tmplWithoutCont->ConvertToLogScale();

    m_ListWithoutCont[r->GetCategory()].push_back( tmplWithoutCont );

    return true;
}

/**
 * Adds the templates in the input path, under the input category.
 */
Bool CTemplateCatalog::Add( const char* templatePath, const std::string& category )
{
    if ( !exists( templatePath ) )
      {
	Log.LogError ( "Path %s does not exist.", templatePath );
	return false;
      }

    path p( templatePath );
    path name = p.leaf();

    std::shared_ptr<CTemplate> tmpl = std::shared_ptr<CTemplate>( new CTemplate( name.c_str(), category ) );

    CSpectrumIOGenericReader asciiReader;
    if( !asciiReader.Read( templatePath, *tmpl ) ) {
        Log.LogError( "Fail to read template: %s", templatePath );
        return false;
    }

    Add( tmpl );
    return true;
}

/**
 * Loads input directory as a collection of categories of templates.
 */
Bool CTemplateCatalog::Load( const char* dirPath )
{
    path pathFound;

    if ( !exists( dirPath ) )
      {
	Log.LogError ( "Path %s does not exist.", dirPath );
        return false;
      }

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
    return true;
}

/**
 * Load files in dirPath under the input category.
 */
Bool CTemplateCatalog::LoadCategory( const path& dirPath, const std::string&  category )
{
    directory_iterator end_itr;
    for ( directory_iterator itr( dirPath ); itr != end_itr; ++itr )
    {
        if ( !is_directory( itr->status() ) )
        {
            if( !Add( itr->path().c_str(), category ) )
                return false;
        }
    }
    return true;
}

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

IMPLEMENT_MANAGED_OBJECT( CTemplateCatalog )

CTemplateCatalog::CTemplateCatalog()
{

}

CTemplateCatalog::~CTemplateCatalog()
{

}

const CTemplate& CTemplateCatalog::GetTemplate( const std::string& category, UInt32 i ) const
{
    return *m_List.at(category)[i];
}


const CTemplate& CTemplateCatalog::GetTemplateWithoutContinuum( const std::string& category, UInt32 i ) const
{
    return *m_ListWithoutCont.at(category)[i];
}

TTemplateRefList CTemplateCatalog::GetTemplate( const TStringList& categoryList ) const
{
    TTemplateRefList list;

    for( Int32 i=0;i<categoryList.size(); i++ )
    {
        for (Int32 j=0; j<GetTemplateCount( categoryList[i] ); i++ )
        {
            list.push_back( m_List.at( categoryList[i] )[j] );
        }
    }

    return list;
}

TTemplateRefList CTemplateCatalog::GetTemplateWithoutContinuum( const TStringList& categoryList ) const
{
    TTemplateRefList list;

    for( Int32 i=0;i<categoryList.size(); i++ )
    {
        for (Int32 j=0; j<GetTemplateCount( categoryList[i] ); i++ )
        {
            list.push_back( m_ListWithoutCont.at( categoryList[i] )[j] );
        }
    }

    return list;
}

TStringList CTemplateCatalog::GetCategoryList() const
{
    TStringList l;

    for( auto it = m_List.begin(); it != m_List.end(); it++ ) {
        l.push_back( it->first );
    }

    return l;
}

UInt32 CTemplateCatalog::GetTemplateCount( const std::string& category ) const
{
    return m_List.at(category).size();
}

Bool CTemplateCatalog::Add( CTemplate& r )
{
    if( r.GetCategory().empty() )
        return false;

    m_List[r.GetCategory()].push_back( &r );

    // Compute continuum substracted spectrum
    CRef<CTemplate> tmplWithoutCont = new CTemplate( r.GetName().c_str(), r.GetCategory() );

    *tmplWithoutCont = r;

    CContinuumIrregularSamplingMedian continuum;

    tmplWithoutCont->RemoveContinuum( continuum );
    tmplWithoutCont->ConvertToLogScale();

    m_ListWithoutCont[r.GetCategory()].push_back( tmplWithoutCont );

    return true;
}

Bool CTemplateCatalog::Add( const char* templatePath, const std::string& category )
{
    if ( !exists( templatePath ) )
        return false;

    path p( templatePath );
    path name = p.leaf();

    CRef<CTemplate> tmpl = new CTemplate( name.c_str(), category );

    CSpectrumIOGenericReader asciiReader;
    if( !asciiReader.Read( templatePath, *tmpl ) ) {
        Log.LogError("Fail to read template: %s", templatePath);
        return false;
    }

    Add( *tmpl );
    return true;
}

Bool CTemplateCatalog::Load( const char* dirPath )
{
    path pathFound;

    if ( !exists( dirPath ) )
        return false;

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

Bool CTemplateCatalog::LoadCategory( const path& dirPath, const std::string&  category )
{
    directory_iterator end_itr;
    for ( directory_iterator itr( dirPath ); itr != end_itr; ++itr )
    {
        if ( !is_directory( itr->status() ) )
        {
            if( ! Add( itr->path().c_str(), category ) )
                return false;
        }
    }
    return true;
}

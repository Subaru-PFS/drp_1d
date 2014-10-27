#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/spectrum/io/asciireader.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/continuum/median.h>

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

const CTemplate& CTemplateCatalog::GetTemplate( CTemplate::ECategory category, UInt32 i ) const
{
    return *m_List[category][i];
}


const CTemplate& CTemplateCatalog::GetTemplateWithoutContinuum( CTemplate::ECategory category, UInt32 i ) const
{
    return *m_ListWithoutCont[category][i];
}

UInt32 CTemplateCatalog::GetTemplateCount( CTemplate::ECategory category ) const
{
    return m_List[category].size();
}

Bool CTemplateCatalog::Add( CTemplate& r )
{
    if( r.GetCategory() == CTemplate::nCategory_None )
        return false;

    m_List[r.GetCategory()].push_back( &r );

    // Compute continuum substracted spectrum
    CRef<CTemplate> tmplWithoutCont = new CTemplate( r.GetName().c_str(), r.GetCategory() );

    *tmplWithoutCont = r;

    tmplWithoutCont->RemoveContinuum<CContinuumMedian>();
    tmplWithoutCont->ConvertToLogScale();

    m_ListWithoutCont[r.GetCategory()].push_back( tmplWithoutCont );

    return true;
}

Bool CTemplateCatalog::Add( const char* templatePath, CTemplate::ECategory category )
{
    if ( !exists( templatePath ) )
        return false;

    path p( templatePath );
    path name = p.leaf();

    CRef<CTemplate> tmpl = new CTemplate( name.c_str(), category );

    CSpectrumIOAsciiReader asciiReader;
    if( !asciiReader.Read( templatePath, *tmpl ) )
        return false;

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
            CTemplate::ECategory category = ConvertStringToCategory( (*it).generic_string() );
            if( category != CTemplate::nCategory_None )
                LoadCategory( itr->path(), category );
        }
    }
    return true;
}

Bool CTemplateCatalog::LoadCategory( const path& dirPath, CTemplate::ECategory category )
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

CTemplate::ECategory CTemplateCatalog::ConvertStringToCategory( const string& category )
{
    if( category == "star" )
        return CTemplate::nCategory_Star;
    else if( category == "galaxy" )
        return CTemplate::nCategory_Galaxy;
    else if( category == "qso" )
        return CTemplate::nCategory_Qso;
    else if( category == "emission" )
        return CTemplate::nCategory_Emission;

    return CTemplate::nCategory_None;
}


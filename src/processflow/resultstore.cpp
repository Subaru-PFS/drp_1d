#include <epic/redshift/processflow/resultstore.h>

#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/template.h>

#include <fstream>


using namespace NSEpic;

namespace bfs = boost::filesystem;

IMPLEMENT_MANAGED_OBJECT( COperatorResultStore );

COperatorResultStore::COperatorResultStore()
{

}

COperatorResultStore::~COperatorResultStore()
{

}

void COperatorResultStore::StoreResult( TResultsMap& map, const std::string& path, const std::string& name, const COperatorResult& result )
{
    std::string scopedName;
    if( ! path.empty() ) {
        scopedName = path;
        scopedName.append( "." );
    }
    scopedName.append( name );

    TResultsMap::iterator it = map.find( name );
    if( it != map.end() )
    {
        DebugError( "Resulat already exist");
        return;
    }

    map[ scopedName ] = &result;
}

Void COperatorResultStore::StorePerTemplateResult( const CTemplate& t, const std::string& path, const std::string& name, const COperatorResult& result )
{
    TPerTemplateResultsMap::iterator it = m_PerTemplateResults.find( t.GetName() );
    if( it == m_PerTemplateResults.end() )
    {
        m_PerTemplateResults[ t.GetName() ] = TResultsMap();
    }

    StoreResult( m_PerTemplateResults[ t.GetName() ], path, name, result );
}

Void COperatorResultStore::StoreGlobalResult( const std::string& path, const std::string& name, const COperatorResult& result )
{
    StoreResult( m_GlobalResults, path, name, result );
}

const COperatorResult* COperatorResultStore::GetPerTemplateResult( const CTemplate& t, const std::string& name ) const
{
    TPerTemplateResultsMap::const_iterator it1 = m_PerTemplateResults.find( t.GetName() );
    if( it1 != m_PerTemplateResults.end() )
    {
        const TResultsMap& m = (*it1).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            return (*it2).second;
        }
    }

    return NULL;
}

TOperatorResultMap COperatorResultStore::GetPerTemplateResult( const std::string& name ) const
{
    TOperatorResultMap map;
    TPerTemplateResultsMap::const_iterator it;
    for( it = m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); ++it )
    {
        std::string tplName = (*it).first;

        const TResultsMap& m = (*it).second;
        TResultsMap::const_iterator it2 = m.find( name );
        if( it2 != m.end() )
        {
            map[tplName] = (*it2).second;
        }
    }

    return map;
}

const COperatorResult* COperatorResultStore::GetGlobalResult( const std::string& name ) const
{
    TResultsMap::const_iterator it = m_GlobalResults.find( name );
    if( it != m_GlobalResults.end() )
    {
        return (*it).second;
    }

    return NULL;
}

void COperatorResultStore::CreateResultStorage( std::fstream& stream, const bfs::path& path, const bfs::path& baseDir ) const
{
    bfs::path outputFilePath = bfs::path( baseDir ).append( path.string() );
    std::fstream outputFile;

    if( bfs::exists( outputFilePath.parent_path() ) == false )
        bfs::create_directories( outputFilePath.parent_path() );

    stream.open( outputFilePath.string().c_str(), std::fstream::out | std::fstream::app);
    if( stream.rdstate() & std::ios_base::failbit )
    {
        return ;
    }

}

void COperatorResultStore::SaveRedshiftResult( const CDataStore& store, const char* dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        CreateResultStorage( outputStream, bfs::path( "redshift.csv" ), bfs::path( dir ) );

        CConstRef<COperatorResult>  result = GetGlobalResult( "redshiftresult" );
        if(result){
            result->SaveLine( store, outputStream );
        }
    }
}

void COperatorResultStore::SaveRedshiftResultHeader(  const char* dir )
{
    // Append best redshift result line to output file
    {
        std::fstream outputStream;
        // Save result at root of output directory
        CreateResultStorage( outputStream, bfs::path( "redshift.csv" ), bfs::path( dir ) );


        outputStream <<  "#Spectrum\tRedshifts\tMerit\tTemplate\tMethod/Path"<< std::endl;
    }
}

Void COperatorResultStore::SaveAllResults( const CDataStore& store, const char* dir ) const
{
    // Store global result
    {
        TResultsMap::const_iterator it;
        for( it=m_GlobalResults.begin(); it != m_GlobalResults.end(); it++ )
        {
            std::string resultName = (*it).first;
            CConstRef<COperatorResult>  result = (*it).second;

            std::fstream outputStream;
            // Save result at root of output directory
            CreateResultStorage( outputStream, bfs::path( resultName + ".csv"), bfs::path( dir ).append( m_SpectrumName ) );
            result->Save( store, outputStream );
        }
    }

    // Store per template results
    {
        TPerTemplateResultsMap::const_iterator it;
        for( it=m_PerTemplateResults.begin(); it != m_PerTemplateResults.end(); it++ )
        {
            std::string templateName = (*it).first;
            const TResultsMap& resultMap = (*it).second;

            TResultsMap::const_iterator it2;
            for( it2=resultMap.begin(); it2 != resultMap.end(); it2++ )
            {
                std::string resultName = (*it2).first;
                CConstRef<COperatorResult>  result = (*it2).second;

                std::fstream outputStream;
                // Save result in sub directories of output directory
                CreateResultStorage( outputStream, bfs::path( templateName ).append( resultName + ".csv"), bfs::path( dir ).append( m_SpectrumName ) );
                result->Save( store, outputStream );
            }
        }
    }
}


std::string COperatorResultStore::GetScope(CConstRef<COperatorResult>  result) const
{
    std::string n="";

    TResultsMap::const_iterator it;
    for( it=m_GlobalResults.begin(); it != m_GlobalResults.end(); it++ )
    {
        CConstRef<COperatorResult>  r = (*it).second;
        if(result == r){
            std::string s = (*it).first;

            std::size_t found = s.rfind(".");
            if (found!=std::string::npos){
                n = s.substr(0,found);
                n = n.append(".");
            }
            break;
        }
    }

    return n;
}

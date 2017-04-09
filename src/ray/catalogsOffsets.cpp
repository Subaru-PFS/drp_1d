#include <epic/core/log/log.h>
#include <epic/redshift/ray/catalogsOffsets.h>

#include <algorithm>    // std::sort
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <string>
#include <fstream>
#include <iostream>

namespace bfs = boost::filesystem;
using namespace NSEpic;
using namespace std;
using namespace boost;


CLineCatalogsOffsets::CLineCatalogsOffsets()
{
    m_Catalogs_relpath = "linecatalogs_offsets/offsetsCatalogs_20170410_m300"; //path to the fixed offset catalog
    Log.LogInfo( "CLineCatalogsOffsets - directory : %s", m_Catalogs_relpath.c_str());
}

CLineCatalogsOffsets::~CLineCatalogsOffsets()
{

}

Bool CLineCatalogsOffsets::SetCtlgRelPath( const char* relPath )
{
    m_Catalogs_relpath = relPath;
    return true;
}

Bool CLineCatalogsOffsets::Init( std::string calibrationPath)
{
    bfs::path calibrationFolder( calibrationPath.c_str() );
    std::string dirPath = (calibrationFolder/m_Catalogs_relpath.c_str()).string();

    bool ret = Load(dirPath.c_str());
    if(!ret)
    {
        Log.LogError("Unable to load the offset catalogs. aborting...");
        return false;
    }else{
        Log.LogInfo("Loaded %d lines offsets catalogs", m_OffsetsCatalog.size());
    }
    return true;
}

Bool CLineCatalogsOffsets::Load( const char* dirPath )
{
    m_OffsetsCatalog.clear();

    //load the catalogs list from the directory
    namespace fs = boost::filesystem;
    fs::path catalogDir(dirPath);

    fs::directory_iterator end_iter;
    std::vector<std::string> catalogList;
    if ( fs::exists(catalogDir) && fs::is_directory(catalogDir))
    {
      for( fs::directory_iterator dir_iter(catalogDir) ; dir_iter != end_iter ; ++dir_iter)
      {
        if (fs::is_regular_file(dir_iter->status()) )
        {
          catalogList.push_back(dir_iter->path().c_str());
        }
      }
    }
    if(catalogList.size()<1)
    {
        return false;
    }
    Log.LogDebug( "CLineCatalogsOffsets - Found %d offsets catalogs", catalogList.size());

    //Load the catalogs in the list
    for(Int32 k=0; k<catalogList.size(); k++)
    {
        LoadCatalog( catalogList[k].c_str() );
    }

    return true;
}

Bool CLineCatalogsOffsets::LoadCatalog( const char* filePath )
{

    SOffsetsCatalog newCatalog;

    ifstream file;
    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
    {
        return false;
    }

    std::string catalogPath = filePath;
    newCatalog.filePath = catalogPath;

    string line;
    // Read file line by line
    while( getline( file, line ) )
    {
        // remove comments
        if(line.compare(0,1,"#",1)==0){
            continue;
        }
        char_separator<char> sep(" \t");

        // Tokenize each line
        typedef tokenizer< char_separator<char> > ttokenizer;
        ttokenizer tok( line, sep );

        // Check if it's not a comment
        ttokenizer::iterator it = tok.begin();
        if( it != tok.end() && *it != "#" )
        {
            string name;
            if( it != tok.end() )
            {
                name = *it;
            }
            else
            {
                return false;
            }

            Float64 offset;
            ++it;
            if( it != tok.end() )
            {
                try
                {
                    offset = lexical_cast<double>(*it);
                }
                catch (bad_lexical_cast)
                {
                    Log.LogError( "Unable to read offset value from file, aborting" );
                    return false;
                }
            }

            newCatalog.Offsets.push_back(offset);
            newCatalog.Names.push_back(name);

        }
    }
    file.close();

    m_OffsetsCatalog.push_back(newCatalog);

    return true;
}

Bool CLineCatalogsOffsets::SetLinesOffsets(CLineModelElementList &LineModelElementList, Int32 index)
{
    Int32 nLines = m_OffsetsCatalog[index].Offsets.size();
    if(index>=nLines)
    {
        return false;
    }

    //loop the offsets in the catalog
    for(Int32 kL=0; kL<nLines; kL++)
    {
        Float64 offset = m_OffsetsCatalog[index].Offsets[kL];
        std::string name = m_OffsetsCatalog[index].Names[kL];
        //find line in the elementList
        for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
        {
            Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
            for(UInt32 j=0; j<nRays; j++){

                if(LineModelElementList.m_Elements[iElts]->m_Rays[j].GetName() == name)
                {
                    LineModelElementList.m_Elements[iElts]->m_Rays[j].SetOffset(offset);
                }

            }


        }

    }
    return true;
}


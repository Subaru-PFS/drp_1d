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


CRayCatalogsOffsets::CRayCatalogsOffsets()
{
    m_Catalog_relpath = "linecatalogs_offsets/linecatalogoffset_uvstackspfs10k_20170407.txt";
    Log.LogInfo( "CRayCatalogsOffsets - Loaded catalog : %s", m_Catalog_relpath.c_str());
}

CRayCatalogsOffsets::~CRayCatalogsOffsets()
{

}

Bool CRayCatalogsOffsets::SetCtlgRelPath( const char* relPath )
{
    m_Catalog_relpath = relPath;
    return true;
}

Bool CRayCatalogsOffsets::Init( std::string calibrationPath)
{
    bfs::path calibrationFolder( calibrationPath.c_str() );
    std::string dirPath = (calibrationFolder/m_Catalog_relpath.c_str()).string();

    bool ret = Load(dirPath.c_str());
    if(!ret)
    {
        Log.LogError("Unable to load the offset catalog. aborting...");
        return false;
    }else{
        Log.LogInfo("Loaded %d lines offsets", m_Offsets.size());
    }
    return true;
}


Bool CRayCatalogsOffsets::Load( const char* filePath )
{

    // Clear current catalog list
    m_Offsets.clear();
    m_Names.clear();

    ifstream file;
    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
    {
        return false;
    }

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

            m_Offsets.push_back(offset);
            m_Names.push_back(name);

        }
    }
    file.close();

    return true;
}

Bool CRayCatalogsOffsets::SetLinesOffsets(CLineModelElementList &LineModelElementList)
{
    //first set all offsets to 0.0 ?

    //loop the offsets in the catalog
    Int32 nLines = m_Offsets.size();
    for(Int32 kL=0; kL<nLines; kL++)
    {
        Float64 offset = m_Offsets[kL];
        std::string name = m_Names[kL];
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


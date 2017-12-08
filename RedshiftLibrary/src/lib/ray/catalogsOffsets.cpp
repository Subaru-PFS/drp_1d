#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/catalogsOffsets.h>

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
    //m_Catalogs_relpath = "linecatalogs_offsets/offsetsCatalogs_20171128_fit"; //path to the fit offset catalog, used only for Linemeas for now
    //m_Catalogs_relpath = "linecatalogs_offsets/offsetsCatalogs_20170410_0"; //path to the fixed offset catalog
    m_Catalogs_relpath = "linecatalogs_offsets/offsetsCatalogs_20170410_m150"; //path to the fixed offset catalog
    //m_Catalogs_relpath = "linecatalogs_offsets/offsetsCatalogs_20170410_m300"; //path to the fixed offset catalog
    //m_Catalogs_relpath = "linecatalogs_offsets/offsetsCatalogs_20170410_steidel"; //path to the steidel offset as of 2017-02 pfs10k simus.
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
    m_Calibration_path = calibrationPath;
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
                    Log.LogError( "Unable to read offset value from offset file, aborting" );
                    return false;
                }
            }

            std::string fitMode="fixed";
            ++it;
            if( it != tok.end() )
            {
                try
                {
                    fitMode = *it;
                }
                catch (bad_lexical_cast)
                {
                    Log.LogError( "Unable to read fittingmode value from offset file, aborting" );
                    return false;
                }
            }

            newCatalog.FittingMode.push_back(fitMode);
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
    //first reset all offsets
    for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
    {
        Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
        for(UInt32 j=0; j<nRays; j++){
            LineModelElementList.m_Elements[iElts]->m_Rays[j].SetOffset(0.0);
            LineModelElementList.m_Elements[iElts]->m_Rays[j].EnableOffsetFit(false); //default is to not fit the offset
        }
    }


    //loop the offsets in the catalog
    for(Int32 kL=0; kL<nLines; kL++)
    {
        Float64 offset = m_OffsetsCatalog[index].Offsets[kL];
        bool enableOffsetFit = m_OffsetsCatalog[index].FittingMode[kL]=="fit";
        std::string name = m_OffsetsCatalog[index].Names[kL];
        //find line in the elementList
        for( UInt32 iElts=0; iElts<LineModelElementList.m_Elements.size(); iElts++ )
        {
            Int32 nRays = LineModelElementList.m_Elements[iElts]->GetSize();
            for(UInt32 j=0; j<nRays; j++){

                if(LineModelElementList.m_Elements[iElts]->m_Rays[j].GetName() == name)
                {
                    LineModelElementList.m_Elements[iElts]->m_Rays[j].SetOffset(offset);
                    LineModelElementList.m_Elements[iElts]->m_Rays[j].EnableOffsetFit(enableOffsetFit);
                    Log.LogInfo( "Offsets: setting line: %s\t offset:%.1f\t fitMode:%s\t fitEnabled:%d ", name.c_str(), offset, m_OffsetsCatalog[index].FittingMode[kL].c_str(), enableOffsetFit);

                }

            }


        }

    }
    return true;
}

Bool CLineCatalogsOffsets::SetLinesOffsetsAutoSelectStack(CLineModelElementList &LineModelElementList, std::string spectrumName)
{
    Int32 offsetCtlgIndex = AutoSelectStackFromReferenceFile(spectrumName);
    if(offsetCtlgIndex>=0)
    {
        bfs::path ctlgPath( m_OffsetsCatalog[offsetCtlgIndex].filePath.c_str() );
        std::string ctlgFileStackNameWExt = ctlgPath.filename().string();
        Log.LogInfo( "CLineCatalogsOffsets: AutoSetUVStack from = %s", ctlgFileStackNameWExt.c_str() );
        SetLinesOffsets(LineModelElementList, offsetCtlgIndex);
    }else
    {
        Log.LogWarning( "CLineCatalogsOffsets: FAILED to AutoSetUVStack with name = %s", spectrumName.c_str() );
        return false;
    }

    return true;
}


Int32 CLineCatalogsOffsets::AutoSelectStackFromReferenceFile(std::string spectrumName)
{
    std::string stack_name = "";

    std::string reference_stacks_relpath = "linecatalogs_offsets/reference_pfs10k_uvstacks.txt";
    bfs::path calibrationFolder( m_Calibration_path.c_str() );
    std::string filePath = (calibrationFolder/reference_stacks_relpath.c_str()).string();

    ifstream file;
    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit )
    {
        return -1;
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
                return -1;
            }
            //check string inclusion
            std::size_t foundstra = spectrumName.find(name.c_str());
            if (foundstra==std::string::npos){
                continue;
            }

            ++it;
            if( it != tok.end() )
            {
                try
                {
                    stack_name = (*it);
                }
                catch (bad_lexical_cast)
                {
                    Log.LogError( "Unable to read stackname value from file, aborting" );
                    return -1;
                }
            }

            //stack name found, break
            break;

        }
    }
    file.close();

    if(stack_name=="")
    {
        return -1;
    }
    //find the corresponding index
    Int32 nCtlg = m_OffsetsCatalog.size();
    Int32 iCtlg = -1;
    for(Int32 k=0; k<nCtlg; k++)
    {
        bfs::path ctlgPath( m_OffsetsCatalog[k].filePath.c_str() );
        std::string ctlgFileStackNameWExt = ctlgPath.filename().string();
        size_t lastindex = ctlgFileStackNameWExt.find_last_of(".");
        std::string ctlgFileStackName = ctlgFileStackNameWExt.substr(0, lastindex);
        if(stack_name==ctlgFileStackName)
        {
           iCtlg = k;
        }
    }

    return iCtlg;
}



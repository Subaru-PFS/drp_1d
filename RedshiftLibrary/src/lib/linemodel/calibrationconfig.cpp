#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/linemodel/calibrationconfig.h>

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


CCalibrationConfigHelper::CCalibrationConfigHelper()
{
    calibration_config_file_relpath = "calibration-config.txt";
}

CCalibrationConfigHelper::~CCalibrationConfigHelper()
{
}

Bool CCalibrationConfigHelper::Init( std::string calibrationPath)
{
    bfs::path calibrationFolder( calibrationPath.c_str() );
    //std::string dirPath = (calibrationFolder.append( tplshapedcatalog_relpath.c_str() )).string();
    std::string filePath = (calibrationFolder/calibration_config_file_relpath.c_str()).string();

    bool ret = Load(filePath.c_str());
    if(!ret)
    {
        Log.LogError("Unable to load the calibration config. aborting...");
        return false;
    }
    return true;
}


Bool CCalibrationConfigHelper::Load( const char* filePath )
{
    ifstream file;
    file.open( filePath, ifstream::in );
    if( file.rdstate() & ios_base::failbit ){
        return false;
    }
    string line;

    // Read file line by line
    Int32 readNums = 0;
    while( getline( file, line ) )
    {
        std::size_t foundstra;
        std::size_t foundstrEqual;

        std::string tplratio_key = "linemodel-tplratio-dir";
        foundstra = line.find(tplratio_key.c_str());
        foundstrEqual = line.find("=");
        if (foundstra!=std::string::npos){
            m_linemodelTplratio_relpath = line.substr(foundstrEqual+1, line.size());
            readNums++;
        }
        std::string offsets_key = "linemodel-offsets-dir";
        foundstra = line.find(offsets_key.c_str());
        foundstrEqual = line.find("=");
        if (foundstra!=std::string::npos){
            m_linemodelOffset_relpath = line.substr(foundstrEqual+1, line.size());
            readNums++;
        }


    }
    file.close();
    if(readNums!=2) //reading 1. tplratiodir, 2. offsetsdir
    {
        return false;
    }

    return true;
}


std::string CCalibrationConfigHelper::Get_linemodelTplratio_relpath()
{
    return m_linemodelTplratio_relpath;
}

std::string CCalibrationConfigHelper::Get_linemodelOffset_relpath()
{
    return m_linemodelOffset_relpath;
}

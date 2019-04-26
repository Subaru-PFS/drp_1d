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

void CCalibrationConfigHelper::Init( std::string calibrationPath)
{
    bfs::path calibrationFolder( calibrationPath.c_str() );
    //std::string dirPath = (calibrationFolder.append( tplshapedcatalog_relpath.c_str() )).string();
    std::string filePath = (calibrationFolder/calibration_config_file_relpath.c_str()).string();

    Load(filePath.c_str());
}


void CCalibrationConfigHelper::Load( const char* filePath )
{
    std::ifstream file;
    file.open( filePath, std::ifstream::in );
    if( file.rdstate() & ios_base::failbit ){
      Log.LogError("Can't load calibration config file [%s]", filePath);
      throw runtime_error("Can't load calibration config file");
    }
    string line;

    // Read file line by line
    Int32 readNums = 0;
    while( getline( file, line ) )
    {
        std::size_t foundstra;
        std::size_t foundstrEqual;

        std::string star_templates_key = "star-templates-dir";
        foundstra = line.find(star_templates_key.c_str());
        foundstrEqual = line.find("=");
        if (foundstra!=std::string::npos){
            m_starTemplates_relpath = line.substr(foundstrEqual+1, line.size());
            readNums++;
        }
        std::string qso_templates_key = "qso-templates-dir";
        foundstra = line.find(qso_templates_key.c_str());
        foundstrEqual = line.find("=");
        if (foundstra!=std::string::npos){
            m_qsoTemplates_relpath = line.substr(foundstrEqual+1, line.size());
            readNums++;
        }


    }
    file.close();
    if(readNums!=2) //reading 1.starstemplates, 2.qsotemplates, ...
    {
      Log.LogError("Invalid calibration config file [%s]", filePath);
      throw runtime_error("Invalid calibration config file");
    }
}

std::string CCalibrationConfigHelper::Get_starTemplates_relpath()
{
    return m_starTemplates_relpath;
}

std::string CCalibrationConfigHelper::Get_qsoTemplates_relpath()
{
    return m_qsoTemplates_relpath;
}


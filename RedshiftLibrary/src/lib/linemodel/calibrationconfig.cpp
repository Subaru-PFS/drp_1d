// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/linemodel/calibrationconfig.h"

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


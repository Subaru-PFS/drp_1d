#ifndef _LINEMODEL_CALIBRATIONCONFIG_
#define _LINEMODEL_CALIBRATIONCONFIG_

#include <RedshiftLibrary/common/datatypes.h>

#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * /ingroup Redshift

 */
class CCalibrationConfigHelper
{

public:
    CCalibrationConfigHelper();
    ~CCalibrationConfigHelper();

    Void Init(std::string calibrationPath);
    Void Load(const char* filePath);

    std::string Get_linemodelTplratio_relpath();
    std::string Get_linemodelOffset_relpath();
    std::string Get_starTemplates_relpath();


private:

    std::string calibration_config_file_relpath="";

    std::string m_linemodelTplratio_relpath="";
    std::string m_linemodelOffset_relpath="";
    std::string m_starTemplates_relpath="";

};


}

#endif

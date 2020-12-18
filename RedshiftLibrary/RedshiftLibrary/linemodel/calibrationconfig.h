#ifndef _REDSHIFT_LINEMODEL_CALIBRATIONCONFIGHELPER_
#define _REDSHIFT_LINEMODEL_CALIBRATIONCONFIGHELPER_

#include <RedshiftLibrary/common/datatypes.h>

#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * \ingroup Redshift
 */
class CCalibrationConfigHelper
{

public:
    CCalibrationConfigHelper();
    ~CCalibrationConfigHelper();

    void Init(std::string calibrationPath);
    void Load(const char* filePath);

    std::string Get_starTemplates_relpath();
    std::string Get_qsoTemplates_relpath();


private:

    std::string calibration_config_file_relpath="";
    std::string m_starTemplates_relpath="";
    std::string m_qsoTemplates_relpath="";

};


}

#endif

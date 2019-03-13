#ifndef _REDSHIFT_PRIORHELPER_CONTINUUM_
#define _REDSHIFT_PRIORHELPER_CONTINUUM_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * /ingroup Redshift

 */
class CPriorHelperContinuum
{

public:

    struct SPriorTZE{
        Float64 logpriorTZE;
        Float64 A_sigma;
        Float64 A_mean;
    };
    typedef std::vector<std::vector<SPriorTZE>> TPriorZEList;
    typedef std::vector<TPriorZEList> TPriorTZEList;

    CPriorHelperContinuum();
    ~CPriorHelperContinuum();

    bool Init( std::string priorDirPath );
    bool LoadFileEZ(const char* filePath , std::vector<std::vector<Float64>>& data);

    bool SetSize(UInt32 size);
    bool SetTNameData(UInt32 k, std::string tname);
    bool SetEZTData(UInt32 k, std::vector<std::vector<Float64>> ezt_data);
    bool SetAGaussmeanData(UInt32 k, std::vector<std::vector<Float64>> agaussmean_data);
    bool SetAGausssigmaData(UInt32 k, std::vector<std::vector<Float64>> agausssigma_data);

    bool GetTplPriorData(std::string tplname, TPriorZEList &zePriorData);

    bool mInitFailed = false;

private:

    TPriorTZEList m_data;
    std::vector<std::string> m_tplnames;

    UInt32 m_nZ = 24;
    UInt32 m_dz = 0.25;
    UInt32 m_z0 = 0.0;

    UInt32 m_nEbv = 10;
    UInt32 m_ebv0 = 0.0;
};

}

#endif

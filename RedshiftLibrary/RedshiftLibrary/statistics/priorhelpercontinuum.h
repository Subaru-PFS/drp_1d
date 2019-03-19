#ifndef _REDSHIFT_PRIORHELPER_CONTINUUM_
#define _REDSHIFT_PRIORHELPER_CONTINUUM_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/range.h>

#include <boost/format.hpp>
#include <float.h>
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
        Float64 priorTZE; //p(ebv, tpl | z)
        Float64 A_sigma;
        Float64 A_mean;
        Float64 betaA;
        Float64 betaTE;
        Float64 betaZ;

        Float64 logprior_precompA; //log(1/sqrt(2 pi sigma_a^2))
        Float64 logprior_precompTE; //log(p(ebv, tpl|z))
        Float64 logprior_precompZ; //log(p(z)/dz)
    };
    typedef std::vector<SPriorTZE> TPriorEList;
    typedef std::vector<TPriorEList> TPriorZEList;
    typedef std::vector<TPriorZEList> TPriorTZEList;

    CPriorHelperContinuum();
    ~CPriorHelperContinuum();

    bool Init( std::string priorDirPath );
    bool LoadFileEZ(const char* filePath , std::vector<std::vector<Float64>>& data);
    bool LoadFileZ(const char* filePath , std::vector<Float64>& data);

    bool SetSize(UInt32 size);
    bool SetTNameData(UInt32 k, std::string tname);
    bool SetEZTData(UInt32 k, std::vector<std::vector<Float64>> ezt_data);
    bool SetAGaussmeanData(UInt32 k, std::vector<std::vector<Float64>> agaussmean_data);
    bool SetAGausssigmaData(UInt32 k, std::vector<std::vector<Float64>> agausssigma_data);
    bool SetPzData(std::vector<Float64> z_data);

    bool GetTplPriorData(std::string tplname,
                         std::vector<Float64> redshifts,
                         TPriorZEList &zePriorData,
                         Int32 outsideZRangeExtensionMode=0);

    bool SetBetaA(Float64 beta);
    bool SetBetaTE(Float64 beta);
    bool SetBetaZ(Float64 beta);

    bool mInitFailed = false;

private:

    TPriorTZEList m_data;
    std::vector<Float64> m_data_pz;
    std::vector<std::string> m_tplnames;

    UInt32 m_nZ = 24;
    Float64 m_dz = 0.25;
    Float64 m_z0 = 0.0;

    UInt32 m_nEbv = 10;
    Float64 m_ebv0 = 0.0;

    Float64 m_betaTE = -1;
    Float64 m_betaA = -1;
    Float64 m_betaZ = -1;

    Float64 m_priorminval = 0.0;//DBL_MIN;
};

}

#endif

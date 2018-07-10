#ifndef _REDSHIFT_LINEMODEL_TEMPLATES_FIT_STORE_
#define _REDSHIFT_LINEMODEL_TEMPLATES_FIT_STORE_


#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>


#include <boost/filesystem.hpp>
#include <vector>

namespace NSEpic
{

class CTemplatesFitStore
{
public:

    struct SValues{
        Float64 redshift;
        Float64 merit;
        Float64 fitAmplitude;
        Float64 fitDustCoeff;
        Int32 fitMeiksinIdx;
        Float64 fitDtM;
        Float64 fitMtM;
        std::string tplName;
    };
    typedef SValues TemplateFitValues;

    CTemplatesFitStore(Float64 minRedshift, Float64 maxRedshift, Float64 stepRedshift, std::string opt_sampling);
    ~CTemplatesFitStore();
    bool Add( Float64 redshift, Float64 merit, Float64 fitAmplitude, Float64 fitDustCoeff, Float64 fitMeiksinIdx, Float64 fitDtM, Float64 fitMtM, std::string tplName );
    std::vector<Float64> GetRedshiftList();
    TemplateFitValues GetFitValues(Float64 redshiftVal);

private:
    std::vector<SValues>    m_fitValues;
    Float64    m_minRedshift;
    Float64    m_maxRedshift;
    Float64    m_stepRedshift;
    std::string m_samplingRedshift;
};



}

#endif

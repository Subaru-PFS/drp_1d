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
        std::string tplName;
        Float64 ismDustCoeff;
        Int32 igmMeiksinIdx;

        Float64 merit;
        Float64 fitAmplitude;
        Float64 fitDtM;
        Float64 fitMtM;
    };
    typedef SValues TemplateFitValues;

    CTemplatesFitStore(Float64 minRedshift, Float64 maxRedshift, Float64 stepRedshift, std::string opt_sampling);
    ~CTemplatesFitStore();

    bool Add(std::string tplName,
             Float64 ismDustCoeff,
             Int32 igmMeiksinIdx,
             Float64 redshift,
             Float64 merit,
             Float64 fitAmplitude,
             Float64 fitDtM,
             Float64 fitMtM);

    void prepareRedshiftList();
    void initFitValues();
    Int32 GetRedshiftIndex(Float64 z);

    std::vector<Float64> GetRedshiftList();
    TemplateFitValues GetFitValues(Float64 redshiftVal, Int32 continuumCandidateRank);

private:
    std::vector<std::vector<SValues>>    m_fitValues; //[nz][n_continuum_candidates]
    Int32 n_continuum_candidates=3;

    Float64    m_minRedshift;
    Float64    m_maxRedshift;
    Float64    m_stepRedshift;
    std::string m_samplingRedshift;

    std::vector<Float64> redshiftgrid;
    std::map<UInt32,UInt32> redshiftgridmap;
    Float64 redshiftgridmapPrecision = 1e-8;
};



}

#endif

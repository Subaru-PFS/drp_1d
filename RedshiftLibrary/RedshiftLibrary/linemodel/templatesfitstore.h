#ifndef _REDSHIFT_LINEMODEL_TEMPLATESFITSTORE_
#define _REDSHIFT_LINEMODEL_TEMPLATESFITSTORE_


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
        Float64 ismEbmvCoeff;
        Int32 igmMeiksinIdx;

        Float64 merit;
        Float64 fitAmplitude;
        Float64 fitAmplitudeError;
        Float64 fitAmplitudeSigma;
        Float64 fitDtM;
        Float64 fitMtM;
        Float64 logprior;
    };
    typedef SValues TemplateFitValues;

    CTemplatesFitStore(const TFloat64List& redshifts);
    ~CTemplatesFitStore();

    bool Add(std::string tplName,
             Float64 ismEbmvCoeff,
             Int32 igmMeiksinIdx,
             Float64 redshift,
             Float64 merit,
             Float64 fitAmplitude,
             Float64 fitAmplitudeError,
             Float64 fitAmplitudeSigma,
             Float64 fitDtM,
             Float64 fitMtM,
             Float64 logprior);

    void initFitValues();
    Int32 GetRedshiftIndex(Float64 z) const;

    const std::vector<Float64> & GetRedshiftList() const;
    TemplateFitValues GetFitValues(Int32 idxz, Int32 continuumCandidateRank) const;
    TemplateFitValues GetFitValues(Float64 redshiftVal, Int32 continuumCandidateRank) const;
    Int32 GetContinuumCount() const;
    Float64 FindMaxAmplitudeSigma(Float64 & z, TemplateFitValues & fitValues);
    //put as public on purpose to avoid the 'old-school' use of getters
    Float64 m_fitContinuum_tplFitSNRMax = 0.0;
    Float64 m_fitContinuum_fitAmplitudeSigmaMAX = 0.0;
    Float64 m_opt_fitcontinuum_maxCount = 2;//default value to 2
private:
    std::vector<std::vector<SValues>>    m_fitValues; //[nz][n_continuum_candidates]
    Int32 n_max_continuum_candidates=10000;
    Int32 n_continuum_candidates=0;

    std::vector<Float64> redshiftgrid;
};



}

#endif

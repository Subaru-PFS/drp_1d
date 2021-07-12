#include "RedshiftLibrary/linemodel/templatesfitstore.h"
#include "RedshiftLibrary/linemodel/elementlist.h"

#include <float.h>

using namespace NSEpic;

CTemplatesFitStore::CTemplatesFitStore(const TFloat64List& redshifts):
    redshiftgrid(redshifts)
{
    initFitValues();
}
/**
 * \brief Empty destructor.
 **/
CTemplatesFitStore::~CTemplatesFitStore()
{
}


/**
 * @brief CTemplatesFitStore::initFitValues
 * allocate the [nz][n continuum candidates values] structure
 */
void CTemplatesFitStore::initFitValues()
{
    for(Int32 kz=0; kz<redshiftgrid.size(); kz++)
    {
        SValues values_unused;
        values_unused.merit = DBL_MAX;
        std::vector<SValues> zfitvals;//(n_max_continuum_candidates, values_unused);
        m_fitValues.push_back(zfitvals);
    }
}

const std::vector<Float64> & CTemplatesFitStore::GetRedshiftList() const 
{
    return redshiftgrid;
}

Int32 CTemplatesFitStore::GetRedshiftIndex(Float64 z) const
{
    auto it = std::find(redshiftgrid.begin(), redshiftgrid.end(), z);
    if (it != redshiftgrid.end())
        return std::distance(redshiftgrid.begin(), it);
    else 
        return -1;
}

/**
 * @brief CTemplatesFitStore::Add
 * @param tplName
 * @param ismEbmvCoeff
 * @param igmMeiksinIdx
 * @param redshift
 * @param merit
 * @param fitAmplitude
 * @param fitAmplitudeError
 * @param fitDtM
 * @param fitMtM
 *
 * brief: try to insert the fit values into the mFitValues table at the correct idxz:
 *   - no insertion if redshift can't be found in the redshiftgrid
 *   - no insertion if the merit is higher than the highest rank continuum candidate
 *   - insertion is done at a given continuum_candidate_rank position wrt merit value
 * @return False if there was a problem.
 */
bool CTemplatesFitStore::Add(std::string tplName,
                             Float64 ismEbmvCoeff,
                             Int32 igmMeiksinIdx,
                             Float64 redshift,
                             Float64 merit,
                             Float64 fitAmplitude,
                             Float64 fitAmplitudeError,
                             Float64 fitAmplitudeSigma,
                             Float64 fitDtM,
                             Float64 fitMtM,
                             Float64 logprior)
{
    SValues tmpSValues;
    tmpSValues.merit = merit;
    tmpSValues.fitAmplitude = fitAmplitude;
    tmpSValues.fitAmplitudeError = fitAmplitudeError;
    tmpSValues.fitAmplitudeSigma = fitAmplitudeSigma;
    tmpSValues.fitDtM = fitDtM;
    tmpSValues.fitMtM = fitMtM;
    tmpSValues.logprior = logprior;
    tmpSValues.ismEbmvCoeff = ismEbmvCoeff;
    tmpSValues.igmMeiksinIdx = igmMeiksinIdx;
    tmpSValues.tplName = tplName;

    //
    Int32 idxz=GetRedshiftIndex(redshift);
    if(idxz<0)
    {
        Log.LogDebug("CTemplatesFitStore::Unable to find z index for redshift=%f",
                     redshift);
        return false;
    }
    //

    //if chi2 val is the lowest, and condition on tplName, insert at position ipos
    Int32 ipos=-1;
    for(Int32 kpos=0; kpos<m_fitValues[idxz].size(); kpos++)
    {
        if(tmpSValues.merit < m_fitValues[idxz][kpos].merit)
        {
            ipos = kpos;
            break;
        }
    }

    if(ipos<0 && m_fitValues[idxz].size()<n_max_continuum_candidates)
    {
        Log.LogDebug("CTemplatesFitStore::Add iz=%d (z=%f) - adding at end position %d (merit=%e, ebmv=%e, imeiksin=%d)",
                     idxz,
                     redshift,
                     m_fitValues[idxz].size(),
                     tmpSValues.merit,
                     tmpSValues.ismEbmvCoeff,
                     tmpSValues.igmMeiksinIdx);

        m_fitValues[idxz].push_back(tmpSValues);

    }else if(ipos<0)
    {
        //nothing to do, merit doesn't qualify the fit result to be stored
    }else{
        Log.LogDebug("CTemplatesFitStore::Add iz=%d (z=%f) - adding at pos=%d (merit=%e, ebmv=%e, imeiksin=%d)",
                     idxz,
                     redshift,
                     ipos,
                     tmpSValues.merit,
                     tmpSValues.ismEbmvCoeff,
                     tmpSValues.igmMeiksinIdx);

        //insert the new SValue and move all the older candidates position according to ipos found
        std::vector<SValues>  tmpBufferValues;
        for(UInt32 ktmp=0; ktmp<m_fitValues[idxz].size(); ktmp++)
        {
            tmpBufferValues.push_back(m_fitValues[idxz][ktmp]);
        }
        SValues values_unused;
        values_unused.merit = DBL_MAX;
        tmpBufferValues.push_back(values_unused);

        UInt32 iOld=0;
        for(UInt32 ktmp=0; ktmp<m_fitValues[idxz].size(); ktmp++)
        {
            if(ipos==ktmp)
            {
                m_fitValues[idxz][ktmp] = tmpSValues;
            }else{
                m_fitValues[idxz][ktmp] = tmpBufferValues[iOld];
                iOld++;
            }
        }
        m_fitValues[idxz].push_back(tmpBufferValues[iOld]);
    }

    //this is not very secure. it should be checked that all redshifts have the same fitValues count
    if(n_continuum_candidates<m_fitValues[idxz].size())
    {
        n_continuum_candidates=m_fitValues[idxz].size();
        Log.LogDebug("CTemplatesFitStore::n_continuum_candidates set to %d)", n_continuum_candidates);
    }


    return true;
}

Int32 CTemplatesFitStore::GetContinuumCount() const
{
    return n_continuum_candidates;
}


CTemplatesFitStore::TemplateFitValues  CTemplatesFitStore::GetFitValues(Int32 idxz, Int32 continuumCandidateRank) const
{
    if(continuumCandidateRank>n_continuum_candidates-1)
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum: candidateRank (%d) >= n_continuum_candidates (%d)",
                     continuumCandidateRank,
                     n_continuum_candidates);
        throw runtime_error("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum");
    }else if(continuumCandidateRank<0)
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum: candidateRank (%d) <0",
                     continuumCandidateRank);
        throw runtime_error("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum");
    }

    if ( (idxz<0) || (idxz > redshiftgrid.size()-1))
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - redshift idx %d is outside range", idxz);        
        throw runtime_error("CTemplatesFitStore::GetFitValues - redshift idx is outside range");
    }
 
    return m_fitValues[idxz][continuumCandidateRank];
}


CTemplatesFitStore::TemplateFitValues  CTemplatesFitStore::GetFitValues(Float64 redshiftVal, Int32 continuumCandidateRank) const
{
    if(continuumCandidateRank>n_continuum_candidates-1)
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum: candidateRank (%d) >= n_continuum_candidates (%d)",
                     continuumCandidateRank,
                     n_continuum_candidates);
        throw runtime_error("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum");

    }else if(continuumCandidateRank<0)
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum: candidateRank (%d) <0",
                     continuumCandidateRank);
        throw runtime_error("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum");
    }

    if(redshiftVal<redshiftgrid[0])
    {
        Log.LogError("CTemplatesFitStore - GetFitValues, looking for redshiftVal=%f, but lt redshiftgrid[0]=%f",
                     redshiftVal,
                     redshiftgrid[0]);
        throw runtime_error("CTemplatesFitStore::GetFitValues - looking for outside range redshiftVal");
    }
    if(redshiftVal>redshiftgrid[redshiftgrid.size()-1])
    {
        Log.LogError("CTemplatesFitStore - GetFitValues, looking for redshiftVal=%f, but ht redshiftgrid[redshiftgrid.size()-1]=%f",
                     redshiftVal,
                     redshiftgrid[redshiftgrid.size()-1]);
        throw runtime_error("CTemplatesFitStore::GetFitValues - looking for outside range redshiftVal");
    }


    //find the idxz
    Int32 idxz=-1;
    idxz = GetRedshiftIndex(redshiftVal);

    if(idxz<0)
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum.");
        Log.LogError("CTemplatesFitStore::GetFitValues - redshiftVal=%e, but lt m_fitValues[0].redshift=%e",
                     redshiftVal,
                     redshiftgrid[0]);
        Log.LogError("CTemplatesFitStore::GetFitValues - redshiftVal=%e, but ht m_fitValues[m_fitValues.size()-1].redshift=%e",
                     redshiftVal,
                     redshiftgrid[redshiftgrid.size()-1]);

        for(Int32 k=0; k<10; k++)
        {
            Log.LogDebug("CTemplatesFitStore::GetFitValues - redshiftVal=%e, lt m_fitValues[%d].redshift=%e",
                         redshiftVal,
                         k,
                         redshiftgrid[k]);
        }
        
        throw runtime_error("CTemplatesFitStore::GetFitValues - cannot find redshiftVal");
    }
    
    return m_fitValues[idxz][continuumCandidateRank];
}

Float64 CTemplatesFitStore::FindMaxAmplitudeSigma(Float64 & z, TemplateFitValues & fitValues)
{
    Int32 icontinuum = 0;
    m_fitContinuum_fitAmplitudeSigmaMAX = -INFINITY;
    //TemplateFitValues fitValues;
    for (Int32 i = 0; i < redshiftgrid.size(); i++)
    {
        const TemplateFitValues &  thisfitValues = m_fitValues[i][icontinuum];
        if (thisfitValues.fitAmplitudeSigma > m_fitContinuum_fitAmplitudeSigmaMAX){
            m_fitContinuum_fitAmplitudeSigmaMAX = thisfitValues.fitAmplitudeSigma;
            z = redshiftgrid[i];
            fitValues = thisfitValues;
        }
    }

    return m_fitContinuum_fitAmplitudeSigmaMAX;

}

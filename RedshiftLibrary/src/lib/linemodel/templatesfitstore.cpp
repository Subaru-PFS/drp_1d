#include <RedshiftLibrary/linemodel/templatesfitstore.h>
#include <RedshiftLibrary/linemodel/elementlist.h>

#include <float.h>

using namespace NSEpic;


CTemplatesFitStore::CTemplatesFitStore(Float64 minRedshift, Float64 maxRedshift, Float64 stepRedshift, std::string opt_sampling)
{
    Float64 marginRedshiftSteps = 3.0;
    m_minRedshift = minRedshift - marginRedshiftSteps*stepRedshift;
    m_maxRedshift = maxRedshift + marginRedshiftSteps*stepRedshift;
    m_stepRedshift = stepRedshift;
    m_samplingRedshift = opt_sampling;

    prepareRedshiftList();

    initFitValues();
}

/**
 * \brief Empty destructor.
 **/
CTemplatesFitStore::~CTemplatesFitStore()
{
}

/**
 * @brief CTemplatesFitStore::prepareRedshiftList
 * prepare the redshift grid and the redshift map
 * @return
 */
void CTemplatesFitStore::prepareRedshiftList()
{
    if(redshiftgridmapPrecision<0.)
    {
        Log.LogError("Failed to Initialize the template fit store redshift map (redshiftgridmapPrecision=%f). Aborting", redshiftgridmapPrecision);
        throw std::runtime_error("Failed to Initialize the template fit store redshift map due to bad redshiftgridmapPrecision value");
    }

    std::vector<Float64> redshifts;
    TFloat64Range redshiftRange = TFloat64Range( m_minRedshift, m_maxRedshift );
    if(m_samplingRedshift=="log")
    {
        redshifts = redshiftRange.SpreadOverLog( m_stepRedshift );
    }else
    {
        redshifts = redshiftRange.SpreadOver( m_stepRedshift );
    }
    redshiftgrid = redshifts;

    for(UInt32 kz=0; kz<redshiftgrid.size(); kz++)
    {
        UInt32 redshiftscaledInt = (UInt32)(redshiftgrid[kz]/redshiftgridmapPrecision+0.5);
        redshiftgridmap.insert(std::make_pair(redshiftscaledInt, kz));
    }
    if(redshiftgridmap.size()!=redshiftgrid.size())
    {
        Log.LogError("Failed to Initialize the template fit store redshift map (n-zgrid=%d, n-zmap=%d). Aborting", redshiftgrid.size(), redshiftgridmap.size());
        throw std::runtime_error("Failed to Initialize the template fit store redshift map");
    }
    return;
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
        std::vector<SValues> zfitvals(n_continuum_candidates, values_unused);
        m_fitValues.push_back(zfitvals);
    }
}

std::vector<Float64> CTemplatesFitStore::GetRedshiftList()
{
    return redshiftgrid;
}

Int32 CTemplatesFitStore::GetRedshiftIndex(Float64 z)
{
    UInt32 redshiftscaledInt = (UInt32)(z/redshiftgridmapPrecision+0.5);
    std::map<UInt32,UInt32>::iterator it = redshiftgridmap.find(redshiftscaledInt);
    if(it != redshiftgridmap.end())
    {
        Int32 outputIdx = (Int32)(it->second);
        return outputIdx;
    }else{
        return -1;
    }
}

/**
 * @brief CTemplatesFitStore::Add
 * @param tplName
 * @param ismDustCoeff
 * @param igmMeiksinIdx
 * @param redshift
 * @param merit
 * @param fitAmplitude
 * @param fitDtM
 * @param fitMtM
 *
 * brief: try to insert the fit values into the mFitValues table at the correct idxz:
 *   - no insertion if redshift can't be found in the redshiftgrid
 *   - no insertion if the merti is higher than the highest rank continuum candidate
 *   - insertion is done at a given continuum_candidate_rank position wrt merit value
 * @return
 */
bool CTemplatesFitStore::Add(std::string tplName,
                             Float64 ismDustCoeff,
                             Int32 igmMeiksinIdx,
                             Float64 redshift,
                             Float64 merit,
                             Float64 fitAmplitude,
                             Float64 fitDtM,
                             Float64 fitMtM)
{
    SValues tmpSValues;
    tmpSValues.merit = merit;
    tmpSValues.fitAmplitude = fitAmplitude;
    tmpSValues.fitDtM = fitDtM;
    tmpSValues.fitMtM = fitMtM;
    tmpSValues.ismDustCoeff = ismDustCoeff;
    tmpSValues.igmMeiksinIdx = igmMeiksinIdx;
    tmpSValues.tplName = tplName;

    //
    Int32 idxz=GetRedshiftIndex(redshift);
    if(idxz<0)
    {
        return false;
    }
    //

    //if chi2 val is the lowest, and condition on tplName, insert at position ipos
    Int32 ipos=-1;
    for(Int32 kpos=0; kpos<n_continuum_candidates; kpos++)
    {
        if(tmpSValues.merit < m_fitValues[idxz][kpos].merit)
        {
            ipos = kpos;
            break;
        }
    }

    if(ipos<0)
    {
        //nothing to do, merit doesn't qualify the fit result to be stored
        return false;
    }else{
        //insert the new SValue and move all the older candidates position according to ipos found
        std::vector<SValues>  tmpBufferValues;
        for(UInt32 ktmp=0; ktmp<m_fitValues[idxz].size(); ktmp++)
        {
            tmpBufferValues.push_back(m_fitValues[idxz][ktmp]);
        }

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
    }
    //

    return true;
}

CTemplatesFitStore::TemplateFitValues CTemplatesFitStore::GetFitValues(Float64 redshiftVal, Int32 continuumCandidateRank)
{
    if(continuumCandidateRank>n_continuum_candidates-1)
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum: candidateRank (%d) >= n_continuum_candidates (%d)",
                     continuumCandidateRank,
                     n_continuum_candidates);
    }else if(continuumCandidateRank<0)
    {
        Log.LogError("CTemplatesFitStore::GetFitValues - cannot find the correct pre-computed continuum: candidateRank (%d) <0",
                     continuumCandidateRank);
    }

    if(redshiftVal<redshiftgrid[0])
    {
        Log.LogError("CTemplatesFitStore - GetFitValues, looking for redshiftVal=%f, but lt redshiftgrid[0]=%f",
                     redshiftVal,
                     redshiftgrid[0]);
        SValues sval;
        sval.tplName="";
        return sval;
    }
    if(redshiftVal>redshiftgrid[redshiftgrid.size()-1])
    {
        Log.LogError("CTemplatesFitStore - GetFitValues, looking for redshiftVal=%f, but ht redshiftgrid[redshiftgrid.size()-1]=%f",
                     redshiftVal,
                     redshiftgrid[redshiftgrid.size()-1]);
        SValues sval;
        sval.tplName="";
        return sval;
    }


    //find the idxz using the zmap
    Int32 idxz=-1;
    if(m_samplingRedshift=="log")
    {
        for(Int32 k=0; k<redshiftgrid.size()-1; k++)
        {
            if(redshiftVal >= redshiftgrid[k] && redshiftVal <= redshiftgrid[k+1] )
            {
                idxz = k;
                break;
            }
        }
    }else{
        idxz = (redshiftVal-m_minRedshift)/m_stepRedshift;
    }

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
            Log.LogError("CTemplatesFitStore::GetFitValues - redshiftVal=%e, but lt m_fitValues[%d].redshift=%e",
                         redshiftVal,
                         k,
                         redshiftgrid[k]);
        }

    }
    return m_fitValues[idxz][continuumCandidateRank];
}

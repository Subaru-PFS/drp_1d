#include <RedshiftLibrary/linemodel/templatesfitstore.h>
#include <RedshiftLibrary/linemodel/elementlist.h>



using namespace NSEpic;


CTemplatesFitStore::CTemplatesFitStore(Float64 minRedshift, Float64 maxRedshift, Float64 stepRedshift, std::string opt_sampling)
{
    Float64 marginRedshiftSteps = 3.0;
    m_minRedshift = minRedshift - marginRedshiftSteps*stepRedshift;
    m_maxRedshift = maxRedshift + marginRedshiftSteps*stepRedshift;
    m_stepRedshift = stepRedshift;
    m_samplingRedshift = opt_sampling;
}

/**
 * \brief Empty destructor.
 **/
CTemplatesFitStore::~CTemplatesFitStore()
{

}

std::vector<Float64> CTemplatesFitStore::GetRedshiftList()
{
    std::vector<Float64> redshifts;
    TFloat64Range redshiftRange = TFloat64Range( m_minRedshift, m_maxRedshift );
    if(m_samplingRedshift=="log")
    {
        redshifts = redshiftRange.SpreadOverLog( m_stepRedshift );
    }else
    {
        redshifts = redshiftRange.SpreadOver( m_stepRedshift );
    }
    return redshifts;
}

bool CTemplatesFitStore::Add(Float64 redshift, Float64 merit, Float64 fitAmplitude, Float64 fitDustCoeff, Float64 fitMeiksinIdx, Float64 fitDtM, Float64 fitMtM , std::string tplName)
{
    SValues values;
    values.redshift = redshift;
    values.merit = merit;
    values.fitAmplitude = fitAmplitude;
    values.fitDustCoeff = fitDustCoeff;
    values.fitMeiksinIdx = fitMeiksinIdx;
    values.fitDtM = fitDtM;
    values.fitMtM = fitMtM;
    values.tplName = tplName;

    m_fitValues.push_back(values);

    return true;
}

CTemplatesFitStore::TemplateFitValues CTemplatesFitStore::GetFitValues(Float64 redshiftVal)
{
    if(redshiftVal<m_fitValues[0].redshift)
    {
        Log.LogError("CTemplatesFitStore - GetFitValues, looking for redshiftVal=%f, but lt m_fitValues[0].redshift=%f", redshiftVal, m_fitValues[0].redshift);
        SValues sval;
        sval.tplName="";
        return sval;
    }
    if(redshiftVal>m_fitValues[m_fitValues.size()-1].redshift)
    {
        Log.LogError("CTemplatesFitStore - GetFitValues, looking for redshiftVal=%f, but ht m_fitValues[m_fitValues.size()-1].redshift=%f", redshiftVal, m_fitValues[m_fitValues.size()-1].redshift);
        SValues sval;
        sval.tplName="";
        return sval;
    }


    Int32 idx=-1;
    if(m_samplingRedshift=="log")
    {
        for(Int32 k=0; k<m_fitValues.size()-1; k++)
        {
            if(redshiftVal >= m_fitValues[k].redshift && redshiftVal <= m_fitValues[k+1].redshift )
            {
                idx = k;
                break;
            }
        }
    }else{
        idx = (redshiftVal-m_minRedshift)/m_stepRedshift;
    }


    return m_fitValues[idx];
}

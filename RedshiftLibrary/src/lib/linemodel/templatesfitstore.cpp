#include <RedshiftLibrary/linemodel/templatesfitstore.h>
#include <RedshiftLibrary/linemodel/elementlist.h>



using namespace NSEpic;


CTemplatesFitStore::CTemplatesFitStore(Float64 minRedshift, Float64 maxRedshift, Float64 stepRedshift)
{
    m_minRedshift = minRedshift;
    m_maxRedshift = maxRedshift;
    m_stepRedshift = stepRedshift;
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
    Int32 nRedshifts = Int32((m_maxRedshift-m_minRedshift)/m_stepRedshift+1);
    for(Int32 k=0; k<nRedshifts; k++)
    {
        Float64 z = m_minRedshift+m_stepRedshift*k;
        redshifts.push_back(z);
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
    Int32 idx = (redshiftVal-m_minRedshift)/m_stepRedshift;

    return m_fitValues[idx];
}

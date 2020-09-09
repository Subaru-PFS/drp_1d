#include <RedshiftLibrary/spectrum/template/template.h>

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <fstream>
#include <iostream>
using namespace NSEpic;
using namespace std;

/**
 * Constructor, empty.
 */
CTemplate::CTemplate( )
{

}

/**
 * Constructor, assigns values to members.
 */
CTemplate::CTemplate( const std::string& name, const std::string& category ) :
  m_Category( category ) 
{
  m_Name = name ;  
}

CTemplate::CTemplate( const std::string& name, const std::string& category,
		      CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis) :
    m_Category( category )
{
  m_Name = name ;
  m_SpectralAxis = spectralAxis;
  m_FluxAxis = fluxAxis;
}

//TODO use copy constructor of CSpectrum
CTemplate::CTemplate( const CTemplate& other):
    m_kDust(other.m_kDust),
    m_meiksinIdx(other.m_meiksinIdx)
{
    m_Name = other.GetName();
    m_FullPath = other.GetFullPath();
    m_SpectralAxis = other.GetSpectralAxis();
    m_FluxAxis = other.GetFluxAxisWithoutIsmIgm();
    if(other.GetFluxAxis().GetSamplesCount())
        m_FluxAxisIsmIgm = other.GetFluxAxis();

    m_estimationMethod = other.GetContinuumEstimationMethod();
    m_dfBinPath = other.GetWaveletsDFBinPath();
    m_medianWindowSize = other.GetMedianWinsize();
    m_nbScales = other.GetDecompScales();

    m_Category = other.GetCategory();
}
//applying the rule of three (Law of the big three in C++11)
CTemplate& CTemplate::operator=(const CTemplate& other)
{
    m_Name = other.GetName();
    m_FullPath = other.GetFullPath();
    m_SpectralAxis = other.GetSpectralAxis();
    m_FluxAxis = other.GetFluxAxisWithoutIsmIgm();
    if(other.GetFluxAxis().GetSamplesCount())
        m_FluxAxisIsmIgm = other.GetFluxAxis();

    m_estimationMethod = other.GetContinuumEstimationMethod();
    m_dfBinPath = other.GetWaveletsDFBinPath();
    m_medianWindowSize = other.GetMedianWinsize();
    m_nbScales = other.GetDecompScales();

    m_Category = other.GetCategory();

    m_kDust = other.m_kDust;
    m_meiksinIdx = other.m_meiksinIdx;
    return *this;
}
/**
 * Destructor, empty.
 */
CTemplate::~CTemplate()
{

}

/**
 * Returns the value stored in m_Name.
 */
const std::string& CTemplate::GetName() const
{
    return m_Name;
}

Int32 CTemplate::GetTemplateByName( const CTemplateCatalog& tplCatalog,
                                    const TStringList& tplCategoryList,
                                    const std::string tplName,
                                    CTemplate& retTpl) 
{
    for( UInt32 i=0; i<tplCategoryList.size(); i++ )
    {
        std::string category = tplCategoryList[i];

        for( UInt32 j=0; j<tplCatalog.GetTemplateCount( category ); j++ )
        {
            const CTemplate& tpl = tplCatalog.GetTemplate( category, j );
            if(tpl.GetName() == tplName){
                retTpl = tpl; //copy constructor
                return 0;
            }
        }
    }
    return -1;
}
/**
 * Returns the value stored in m_Category.
 */
const std::string& CTemplate::GetCategory() const
{
    return m_Category;
}


/*void CTemplate::SetIsmIgmLambdaRange(Int32 kstart, Int32 kend) const
{
    m_kstart = kstart;
    m_kend = kend;
}*/
/**
 * Saves the template in the given filePath.
 */
Bool CTemplate::Save( const char* filePath ) const
{
    std::fstream file;

    file.open( filePath, fstream::out );
    if( file.rdstate() & ios_base::failbit )
    {
        return false;
    }

    const CSpectrumSpectralAxis& spectralAxis = GetSpectralAxis();
    CSpectrumSpectralAxis spectralAxisCopy(spectralAxis);
    bool logScale = spectralAxisCopy.IsInLogScale();
    //alway save in Linear Scale
    if(logScale)
    {
        spectralAxisCopy.ConvertToLinearScale();
    }
    const CSpectrumFluxAxis& fluxAxis = GetFluxAxis();
    for ( Int32 i=0; i<GetSampleCount(); i++)
    {
        file.precision(10);
        file  <<  spectralAxisCopy[i] << "\t" ;

        file.precision(10);
        file<< fluxAxis[i] << std::endl;
    }
    file.close();
    return true;
}

//Calzetti extinction
bool  CTemplate::ApplyDustCoeff(Int32 kDust)
{
    const TAxisSampleList & Xtpl = m_SpectralAxis.GetSamplesVector();

    m_FluxAxisIsmIgm.SetSize(m_SpectralAxis.GetSamplesCount());

    for(Int32 k= 0; k< m_SpectralAxis.GetSamplesCount(); k++)
    {
        if(kDust!= m_kDust){//if we changed the kDust value
            m_kDust = kDust;
            /*m_computedDustCoeff.resize(m_SpectralAxis.GetSamplesCount());
            //reinit to 1
            std::fill(m_computedDustCoeff.begin(), m_computedDustCoeff.end(), 1);*/
            /*if(!m_computedMeiksingCoeff.size()){
                m_computedMeiksingCoeff.resize(m_SpectralAxis.GetSamplesCount());
                std::fill(m_computedMeiksingCoeff.begin(), m_computedMeiksingCoeff.end(), 1);
            }*/
            m_computedDustCoeff[k] = m_ismCorrectionCalzetti.getDustCoeff( kDust, Xtpl[k]);
        } 
        m_FluxAxisIsmIgm[k] = m_FluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
    }
    return true;
}


bool  CTemplate::ApplyMeiksinCoeff(Int32 meiksinIdx, Float64 redshift)
{
    Int32 redshiftIdx = m_igmCorrectionMeiksin.GetRedshiftIndex(redshift); //index for IGM Meiksin redshift range
    Bool igmCorrectionAppliedOnce = false;
    const TAxisSampleList & Xtpl = m_SpectralAxis.GetSamplesVector();

    m_FluxAxisIsmIgm.SetSize(m_SpectralAxis.GetSamplesCount());

    for(Int32 k = 0; k < m_SpectralAxis.GetSamplesCount(); k++){
         
        if(Xtpl[k] <= m_igmCorrectionMeiksin.GetLambdaMax()){
            if(m_meiksinIdx != meiksinIdx){
                m_meiksinIdx = meiksinIdx;
                /*m_computedMeiksingCoeff.resize(m_SpectralAxis.GetSamplesCount());
                std::fill(m_computedMeiksingCoeff.begin(), m_computedMeiksingCoeff.end(), 1);*/
                /*if(!m_computedDustCoeff.size()){
                    m_computedDustCoeff.resize(m_SpectralAxis.GetSamplesCount());
                    std::fill(m_computedDustCoeff.begin(), m_computedDustCoeff.end(), 1);
                }*/
                Int32 kLbdaMeiksin = 0;
                if(Xtpl[k] >= m_igmCorrectionMeiksin.GetLambdaMin())
                {
                    kLbdaMeiksin = Int32(Xtpl[k] - m_igmCorrectionMeiksin.GetLambdaMin());
                }else //if lambda lower than min meiksin value, use lower meiksin value
                {
                    kLbdaMeiksin = 0;
                }
                m_computedMeiksingCoeff[k] = m_igmCorrectionMeiksin.m_corrections[redshiftIdx].fluxcorr[meiksinIdx][kLbdaMeiksin];
            }
            m_FluxAxisIsmIgm[k] = m_FluxAxis[k]*m_computedMeiksingCoeff[k]*m_computedDustCoeff[k];
            igmCorrectionAppliedOnce = true;
        }
    }
    return igmCorrectionAppliedOnce;
}

//reinit ism/igm configuration when we change redshift value
bool CTemplate::ReinitIsmIgmConfig(){
    m_kDust = -1;
    m_meiksinIdx = -1;
    //no more need to copy the fluxaxis since we decided to track the correction application status and to read directly from fluxAxis
    //set size of computedCorreCoeff and fill it with ones
    //m_computedMeiksingCoeff.resize(m_SpectralAxis.GetSamplesCount());
    //std::fill(m_computedMeiksingCoeff.begin(), m_computedMeiksingCoeff.end(), 1);

    //m_computedDustCoeff.resize(m_SpectralAxis.GetSamplesCount());
    //std::fill(m_computedDustCoeff.begin(), m_computedDustCoeff.end(), 1);         
    return true;
}

bool CTemplate::ReinitIsmIgmComputedCoeffs(Int32 kDust, Int32 meiksinIdx){
    if(m_meiksinIdx != meiksinIdx){
        m_computedMeiksingCoeff.resize(m_SpectralAxis.GetSamplesCount());
        std::fill(m_computedMeiksingCoeff.begin(), m_computedMeiksingCoeff.end(), 1);
    }
    if(kDust!= m_kDust){
        m_computedDustCoeff.resize(m_SpectralAxis.GetSamplesCount());
        std::fill(m_computedDustCoeff.begin(), m_computedDustCoeff.end(), 1);
    }
    return true;
}

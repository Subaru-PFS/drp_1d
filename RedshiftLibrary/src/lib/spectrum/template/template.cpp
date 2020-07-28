#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/mask.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <fstream>

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
  this->m_Name = name ;  
}

CTemplate::CTemplate( const std::string& name, const std::string& category,
		      CSpectrumSpectralAxis& spectralAxis, CSpectrumFluxAxis& fluxAxis) :
    m_Category( category )
{
  this->m_Name = name ;
  m_SpectralAxis = spectralAxis;
  m_FluxAxis = fluxAxis;
  m_FluxAxisIsmIgm = fluxAxis;//by default, same as fluxAxis
}

//TODO use copy constructor of CSpectrum
CTemplate::CTemplate( const CTemplate& other)
{
    m_Name = other.GetName();
    m_FullPath = other.GetFullPath();
    m_SpectralAxis = other.GetSpectralAxis();
    m_FluxAxis = other.GetFluxAxis();
    m_FluxAxisIsmIgm = other.GetFluxAxisIsmIgm();

    m_estimationMethod = other.GetContinuumEstimationMethod();
    m_dfBinPath = other.GetWaveletsDFBinPath();
    m_medianWindowSize = other.GetMedianWinsize();
    m_nbScales = other.GetDecompScales();

    m_Category = other.GetCategory();
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


void CTemplate::SetFluxCorrectionIsmIgm(std::shared_ptr<CSpectrumFluxCorrectionCalzetti> ismCorrectionCalzetti, std::shared_ptr<CSpectrumFluxCorrectionMeiksin> igmCorrectionMeiksin) const    
{
    m_ismCorrectionCalzetti = ismCorrectionCalzetti;
    m_igmCorrectionMeiksin = igmCorrectionMeiksin;
}
void CTemplate::SetIsmIgmLambdaRange(Int32 kstart, Int32 kend) const
{
    m_kstart = kstart;
    m_kend = kend;
}
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
bool  CTemplate::ApplyDustCoeff(Int32 kDust) const 
{
    m_kDust = kDust;
    const TAxisSampleList & Xtpl = m_SpectralAxis.GetSamplesVector();
    TAxisSampleList  Ytpl = m_FluxAxisIsmIgm.GetSamplesVector();
    for(Int32 k= m_kstart; k<= m_kend; k++)
    {
        Float64 coeffDust = m_ismCorrectionCalzetti->getDustCoeff( kDust, Xtpl[k]);
        Ytpl[k] *= coeffDust;
    }
    return true;
}


bool  CTemplate::ApplyMeiksinCoeff(Int32 meiksinIdx, Float64 redshift)const 
{
    m_meiksinIdx = meiksinIdx;
    Int32 redshiftIdx = m_igmCorrectionMeiksin->GetRedshiftIndex(redshift); //index for IGM Meiksin redshift range
        
    Float64 coeffIGM = 1.0;
    Bool igmCorrectionAppliedOnce = false;
    const TAxisSampleList & Xtpl = m_SpectralAxis.GetSamplesVector();
    TAxisSampleList  Ytpl = m_FluxAxisIsmIgm.GetSamplesVector();
    for(Int32 k = m_kstart; k <= m_kend; k++){
         
        if(Xtpl[k] <= m_igmCorrectionMeiksin->GetLambdaMax()){
            Int32 kLbdaMeiksin = 0;
            if(Xtpl[k] >= m_igmCorrectionMeiksin->GetLambdaMin())
            {
                kLbdaMeiksin = Int32(Xtpl[k]- m_igmCorrectionMeiksin->GetLambdaMin());
            }else //if lambda lower than min meiksin value, use lower meiksin value
            {
                kLbdaMeiksin = 0;
            }

            coeffIGM = m_igmCorrectionMeiksin->m_corrections[redshiftIdx].fluxcorr[meiksinIdx][kLbdaMeiksin];
            Ytpl[k] *= coeffIGM;
            igmCorrectionAppliedOnce = true;
        }
    }
    return igmCorrectionAppliedOnce;
}

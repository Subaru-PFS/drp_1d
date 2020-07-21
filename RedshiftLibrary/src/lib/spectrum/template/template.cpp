#include <RedshiftLibrary/spectrum/template/template.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/common/mask.h>
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
}

//TODO use copy constructor of CSpectrum
CTemplate::CTemplate( const CTemplate& other)
{
    m_Name = other.GetName();
    m_FullPath = other.GetFullPath();
    m_SpectralAxis = other.GetSpectralAxis();
    m_FluxAxis = other.GetFluxAxis();

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

/**
 * Returns the value stored in m_Category.
 */
const std::string& CTemplate::GetCategory() const
{
    return m_Category;
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

Int32  CTemplate::RebinTemplate( const CSpectrum& spectrum, 
                                Float64 redshift,
                                const TFloat64Range& lambdaRange,
                                std::string opt_interp,
                                //return variables
                                CSpectrumSpectralAxis& spcSpectralAxis_restframe,
                                CSpectrum& itplTplSpectrum,
                                CMask& itplMask,
                                TFloat64Range& currentRange,
                                Float64& overlapRate,
                                Float64 overlapThreshold,
                                TFloat64List& YtplRawBuffer) const
{
    //TFloat64List YtplRawBuffer;
    const CSpectrumSpectralAxis& spcSpectralAxis = spectrum.GetSpectralAxis();
    const CSpectrumFluxAxis& spcFluxAxis = spectrum.GetFluxAxis();

    const CSpectrumSpectralAxis& tplSpectralAxis = this->GetSpectralAxis();
    const CSpectrumFluxAxis& tplFluxAxis = this->GetFluxAxis();

    //CSpectrumSpectralAxis spcSpectralAxis_restframe;

    Float64 onePlusRedshift = 1.0 + redshift;

    //shift lambdaRange backward to be in restframe
    TFloat64Range spcLambdaRange_restframe;
    TFloat64Range lambdaRange_restframe( lambdaRange.GetBegin() / onePlusRedshift,
                                         lambdaRange.GetEnd() / onePlusRedshift );

    //redshift in restframe the tgtSpectralAxis, i.e., division by (1+Z)
    spcSpectralAxis_restframe.ShiftByWaveLength(spcSpectralAxis, onePlusRedshift, CSpectrumSpectralAxis::nShiftBackward);
    spcSpectralAxis_restframe.ClampLambdaRange( lambdaRange_restframe, spcLambdaRange_restframe );
                                         
    // Compute clamped lambda range over template in restframe
    TFloat64Range tplLambdaRange;
    tplSpectralAxis.ClampLambdaRange( lambdaRange_restframe, tplLambdaRange );
    // Compute the intersected range
    TFloat64Range intersectedLambdaRange( 0.0, 0.0 );
    TFloat64Range::Intersect( tplLambdaRange, spcLambdaRange_restframe, intersectedLambdaRange );

    //this refers to the template we are rebinning
    this->Rebin( intersectedLambdaRange, spcSpectralAxis_restframe, itplTplSpectrum, itplMask, opt_interp);   

    CSpectrumFluxAxis& itplTplFluxAxis = itplTplSpectrum.GetFluxAxis();
    const CSpectrumSpectralAxis& itplTplSpectralAxis = itplTplSpectrum.GetSpectralAxis();
    
    //overlapRate
    overlapRate = spcSpectralAxis_restframe.IntersectMaskAndComputeOverlapRate( lambdaRange_restframe, itplMask );

    // Check for overlap rate
    if( overlapRate < overlapThreshold || overlapRate<=0.0 )
    {
        //status = nStatus_NoOverlap; 
        return -1 ;
    }


    const TAxisSampleList & Xtpl = itplTplSpectralAxis.GetSamplesVector();
    TAxisSampleList & Ytpl = itplTplFluxAxis.GetSamplesVector();
    const TAxisSampleList & Xspc = spcSpectralAxis_restframe.GetSamplesVector();
    const TAxisSampleList & Yspc = spcFluxAxis.GetSamplesVector();
    TFloat64Range logIntersectedLambdaRange( log( intersectedLambdaRange.GetBegin() ), log( intersectedLambdaRange.GetEnd() ) );
    //the spectral axis should be in the same scale
    currentRange = logIntersectedLambdaRange;
    if( spcSpectralAxis_restframe.IsInLinearScale() != tplSpectralAxis.IsInLinearScale() )
    {
        Log.LogError( "    chisquare operator: data and model not in the same scale (lin/log) ! Aborting.");
        //status = nStatus_DataError;
        return -2;
    }
    if(spcSpectralAxis_restframe.IsInLinearScale()){
        currentRange = intersectedLambdaRange;
    }

    for(Int32 k=0; k<itplTplSpectralAxis.GetSamplesCount(); k++)
    {
        YtplRawBuffer.push_back(Ytpl[k]);
    }

    return 0;//YtplRawBuffer;
}
//Calzetti extinction
bool  CTemplate::ApplyDustCoeff(Int32 kDust, const TAxisSampleList & Xtpl, TAxisSampleList & Ytpl/*, TFloat64List Ytpl*/, Int32 kstart, Int32 kend, CSpectrumFluxCorrectionCalzetti* ismCorrectionCalzetti) const 
{
    
    for(Int32 k=kstart; k<=kend; k++)
    {
        Float64 coeffDust = ismCorrectionCalzetti->getDustCoeff( kDust, Xtpl[k]);
        Ytpl[k] *= coeffDust;
    }
    return true;
}


bool  CTemplate::ApplyMeiksinCoeff(Int32 meiksinIdx, const TAxisSampleList & Xtpl, TAxisSampleList & Ytpl/*, TFloat64List Ytpl*/, Int32 kstart, Int32 kend, Int32 redshiftIdx, CSpectrumFluxCorrectionMeiksin* igmCorrectionMeiksin)const 
{
                
    Float64 coeffIGM = 1.0;
    Bool igmCorrectionAppliedOnce = false;
    for(Int32 k=kstart; k<=kend; k++){
        if(Xtpl[k] <= igmCorrectionMeiksin->GetLambdaMax()){
            Int32 kLbdaMeiksin = 0;
            if(Xtpl[k] >= igmCorrectionMeiksin->GetLambdaMin())
            {
                kLbdaMeiksin = Int32(Xtpl[k]-igmCorrectionMeiksin->GetLambdaMin());
            }else //if lambda lower than min meiksin value, use lower meiksin value
            {
                kLbdaMeiksin = 0;
            }

            coeffIGM = igmCorrectionMeiksin->m_corrections[redshiftIdx].fluxcorr[meiksinIdx][kLbdaMeiksin];
            Ytpl[k] *= coeffIGM;
            igmCorrectionAppliedOnce = true;
        }
    }
    return igmCorrectionAppliedOnce;
}

#include <epic/redshift/operator/raydetection.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/spectrum/fluxaxis.h>
#include <epic/redshift/spectrum/spectralaxis.h>
#include <epic/redshift/common/median.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/gaussianfit/gaussianfit.h>

#include <math.h>
#include <float.h>
#include <stdio.h>

using namespace NSEpic;
IMPLEMENT_MANAGED_OBJECT(CRayDetection)

CRayDetection::CRayDetection()
{

}

CRayDetection::~CRayDetection()
{

}

const CRayDetectionResult* CRayDetection::Compute( const CSpectrum& spectrum, const TLambdaRange& lambdaRange, const TInt32RangeList& resPeaks, const TInt32RangeList& resPeaksEnlarged )
{
    const CSpectrum& spc = spectrum;
    const CSpectrumFluxAxis fluxAxis = spc.GetFluxAxis();

    // detect possible peaks
    Float64 winsize = 250.0;
    Float64 cut = 5.0;
    Float64 strongcut = 2.0;
    Float64 minsize = 3;
    Float64 maxsize = 90;

    UInt32 nPeaks = resPeaks.size();

    CRayDetectionResult* result = new CRayDetectionResult();

    // filter the peaks with gaussian fit and create the detected rays catalog
    for( UInt32 j=0; j<nPeaks; j++ )
    {
        bool toAdd = true;
        //find gaussian fit
        CGaussianFit fitter;
        CGaussianFit::EStatus status = fitter.Compute( spc, TInt32Range( resPeaksEnlarged[j].GetBegin(), resPeaksEnlarged[j].GetEnd() ) );
        if(status!=NSEpic::CGaussianFit::nStatus_Success){
            continue;
        }

        Float64 gaussAmp;
        Float64 gaussPos;
        Float64 gaussWidth;
        fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
        Float64 gaussCont;
        fitter.GetResultsPolyCoeff0( gaussCont );

        /* check pos
        if(){
            toAdd = false;
        }
        //*/
        // check amp
        if(gaussAmp<0){
            toAdd = false;
        }
        // check width
        if(toAdd){
            Float64 FWHM_FACTOR=2.35;
            if(gaussWidth<0){
                toAdd = false;
            }else{
                Float64 fwhm = FWHM_FACTOR*gaussWidth;
                if(fwhm<minsize){
                    toAdd = false;
                }
                if(fwhm>maxsize){
                    toAdd = false;
                }
            }
        }

        //find max value and pos
        Float64 max_value = DBL_MIN;
        Int32 max_index = -1;
        for( Int32 k=resPeaks[j].GetBegin(); k<resPeaks[j].GetEnd()+1; k++ )
        {
            if(max_value < fluxAxis[k]){
                max_value = fluxAxis[k];
                max_index = k;
            }
        }

        // Check if gaussian fit is very different from peak itself
        if(toAdd){
            // check max_gauss vs max_spectrum
            Float64 gaussAmp_with_cont = gaussAmp + gaussCont;
            if(gaussAmp_with_cont/max_value <= 0.65 || gaussAmp_with_cont/max_value >= 1.35){
                toAdd = false;
            }
            if(fabs(gaussPos-spc.GetSpectralAxis()[max_index])>3.*spc.GetResolution()){
                toAdd = false;
            }

        }

        // check type weak or strong
        Int32 force = 1; //weak by default
        if(toAdd){
            // strong/weak test to do
            const Float64* fluxData = fluxAxis.GetSamples();
            Int32 windowSampleCount = winsize / spc.GetResolution();
            CMedian<Float64> medianProcessor;
            int left = max(0, (Int32)(max_index-windowSampleCount/2.0+0.5) ) ;
            int right = min((Int32)fluxAxis.GetSamplesCount()-1, (Int32)(max_index + windowSampleCount/2.0) )+1;
            int size_odd = right - left +1;
            if( (int)(size_odd/2) == size_odd/2 ){
                size_odd -= 1;
            }

            Float64 med = medianProcessor.Find( fluxData + left, size_odd);
            Float64 xmad = XMadFind( fluxData + left, size_odd , med);
            // TODO: check with the noise spectrum
            /*if noise!=None:
            noise_mean=noise[left_index:right_index].mean()
            if noise_mean>xmadm:
                # Use noise file
                xmadm=noise_mean
                self.msg.debug(self.__module__, '18.96', 0, self.spectrum.name, "%4.1f" % pos)
            */
            Float64 max_value_no_continuum = max_value - med;
            Float64 ratioAmp=max_value_no_continuum/xmad;

            if(ratioAmp<cut){
                toAdd = false;
            }else if(ratioAmp>cut*strongcut){
                force = 2; //strong
            }
        }

        if(toAdd){
            char buffer [64];
            sprintf(buffer,"detected_peak_%d",j);
            std::string peakName = buffer;
            result->RayCatalog.Add( CRay( peakName, gaussPos, CRay::nType_Emission, force ) );

        }
    }

    return result;
}

Float64 CRayDetection::XMadFind( const Float64* x, Int32 n, Float64 median )
{
    std::vector<Float64> xdata;
    Float64 xmadm = 0.0;

    xdata.reserve( n );

    for( Int32 i=0;i<n; i++ )
    {
        xdata[i] = fabs( x[i]-median );
    }

    CQuickSort<Float64> sort;

    sort.Sort( xdata.data(), n);

    if( ((float)n)/2. - int(n/2.) == 0 )
    {
        UInt32 i1 = n/2 -1;
        UInt32 i2 = n/2;
        xmadm = 0.5*(xdata[i1]+xdata[i2]);
    }
    else
    {
        UInt32 i1 = int(n/2);
        xmadm = xdata[i1];
    }

    return xmadm;
}

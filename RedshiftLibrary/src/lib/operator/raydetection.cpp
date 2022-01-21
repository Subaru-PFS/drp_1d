// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include "RedshiftLibrary/debug/assert.h"
#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/common/median.h"
#include "RedshiftLibrary/gaussianfit/gaussianfit.h"
#include "RedshiftLibrary/operator/raydetection.h"
#include "RedshiftLibrary/operator/raydetectionresult.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/fluxaxis.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"

#include "boost/format.hpp"

#include <math.h>
#include <float.h>
#include <stdio.h>

using namespace NSEpic;

/**
 * Attribution constructor.
 */
CLineDetection::CLineDetection( Int32 type, Float64 cut, Float64 strongcut, Float64 winsize, Float64 minsize, Float64 maxsize, bool disableFitQualityCheck )
{
    FWHM_FACTOR = 2.35;

    m_winsize = winsize;
    m_minsize = minsize;
    m_maxsize = maxsize;
    m_cut = cut;
    m_strongcut = strongcut;

    m_disableFitQualityCheck = disableFitQualityCheck;

    m_type = type;

    m_bypassDebug = true;
}

/**
 * Empty destructor.
 */
CLineDetection::~CLineDetection()
{

}

/**
 * Produce a CLineDetectionResult object containing the lines detected in the input peaks lists.
 */
std::shared_ptr<const CLineDetectionResult> CLineDetection::Compute( const CSpectrum& spectrum, const TLambdaRange& lambdaRange, const TInt32RangeList& resPeaks, const TInt32RangeList& resPeaksEnlarged )
{
    const CSpectrum& spc = spectrum;
    const CSpectrumFluxAxis fluxAxis = spc.GetFluxAxis();
    const CSpectrumSpectralAxis spectralAxis = spc.GetSpectralAxis();

    UInt32 nPeaks = resPeaks.size();

    auto result = std::make_shared<CLineDetectionResult>();

    //retest list
    TInt32RangeList retestPeaks;
    TGaussParamsList retestGaussParams;

    // filter the peaks with gaussian fit and create the detected rays catalog
    for( UInt32 j=0; j<nPeaks; j++ )
    {
        bool toAdd = true;
        //find gaussian fit
        CGaussianFit fitter;

        //limit the peakRange
        // optionally limit the size of the fitrange
        //TInt32Range fitRange = LimitGaussianFitStartAndStop( j, resPeaksEnlarged, spectralAxis.GetSamplesCount(), spectralAxis);
        TInt32Range fitRange = resPeaksEnlarged[j];

        CGaussianFit::EStatus status = fitter.Compute( spc, TInt32Range( fitRange.GetBegin(), fitRange.GetEnd() ) );
        if( status!=NSEpic::CGaussianFit::nStatus_Success )
	  {
            std::string status = ( boost::format( "Peak_%1% : Fitting failed" ) % j ).str();
	    if ( ! m_bypassDebug )
	      Log.LogDebug ( "Peak_%d : Fitting failed", j );
            result->PeakListDetectionStatus.push_back( status );
            continue;
	  }

        Float64 gaussAmp;
        Float64 gaussPos;
        Float64 gaussWidth;
        fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
        Float64 gaussAmpErr;
        Float64 gaussPosErr;
        Float64 gaussWidthErr;
        fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );

        Float64 gaussCont;
        fitter.GetResultsPolyCoeff0( gaussCont );

        // check amp
        if( gaussAmp<0 )
	  {
            toAdd = false;
            std::string status = (boost::format( "Peak_%1% : GaussAmp negative" ) % j).str();
	    if ( ! m_bypassDebug )
	      Log.LogDebug ( "Peak_%d : GaussAmp negative", j );
            result->PeakListDetectionStatus.push_back( status );
	  }

        if( toAdd ) // check width
	  {
            if( gaussWidth<0 )
	      {
                toAdd = false;
                std::string status = (boost::format( "Peak_%1% : GaussWidth negative" ) % j).str();
		if ( ! m_bypassDebug )
		  Log.LogDebug ( "Peak_%d : GaussWidth negative", j );
                result->PeakListDetectionStatus.push_back( status );
	      }
	    else
	      {
                Float64 fwhm = FWHM_FACTOR*gaussWidth;
                if( fwhm<m_minsize )
		  {
                    toAdd = false;
                    std::string status = (boost::format( "Peak_%1% : fwhm<m_minsize" ) % j).str();
		    if ( ! m_bypassDebug )
		      Log.LogDebug ( "Peak_%d : fwhm<m_minsize", j );
                    result->PeakListDetectionStatus.push_back( status );
		  }
                if( fwhm>m_maxsize )
		  {
                    toAdd = false;
                    std::string status = (boost::format( "Peak_%1% : fwhm>m_maxsize" ) % j).str();
		    if ( ! m_bypassDebug )
		      Log.LogDebug ( "Peak_%d : fwhm>m_maxsize", j );
                    result->PeakListDetectionStatus.push_back( status );
		  }
	      }
	  } // check width

        if( toAdd && !m_disableFitQualityCheck ) // Check if gaussian fit is very different from peak itself
	  {
            //find max value and pos
            Float64 max_value = -DBL_MAX;
            Int32 max_index = -1;
            for( Int32 k=resPeaks[j].GetBegin(); k<resPeaks[j].GetEnd()+1; k++ )
            {
                if( max_value < fluxAxis[k] )
		  {
                    max_value = fluxAxis[k];
                    max_index = k;
		  }
            }

            // check flux max_gauss vs flux max_raw_spectrum
            Float64 gaussAmp_with_cont = gaussAmp + gaussCont;
            if( gaussAmp_with_cont/max_value<=0.65 || gaussAmp_with_cont/max_value>=1.35 )
	      {
                toAdd = false;
                std::string status = (boost::format( "Peak_%1% : gaussAmp far from spectrum max_value" ) % j).str();
		if ( ! m_bypassDebug )
		  Log.LogDebug ( "Peak_%d : gaussAmp far from spectrum max_value", j );
                result->PeakListDetectionStatus.push_back( status );
	      }

            if( 0 ) // check gaussPos vs position of the max.
	      { // tolerance on the position in terms of number of samples
                //reg. sampling, TAG: IRREGULAR_SAMPLING
                //if(fabs(gaussPos-spc.GetSpectralAxis()[max_index])>3.*spc.GetResolution()){
                //    toAdd = false;
                //}
                //irregular sampling
                Int32 nsamplestol = 3;
                Float64 error3samples = 1.0;
                if( max_index<nsamplestol )
		  {
                    error3samples = spc.GetSpectralAxis()[max_index+nsamplestol];
		  }
		else 
		  {
		    if( max_index>spc.GetSampleCount()-nsamplestol-1 )
		      {
			error3samples = spc.GetSpectralAxis()[max_index-nsamplestol];
		      }
		    else
		      {
			error3samples = (spc.GetSpectralAxis()[max_index+nsamplestol]-spc.GetSpectralAxis()[max_index-nsamplestol])/2.0;
		      }
		  }
                Float64 diffPos = fabs(gaussPos-spc.GetSpectralAxis()[max_index]);
                if( diffPos > error3samples )
		  {
                    toAdd = false;
                    std::string status = (boost::format( "Peak_%1% : gaussPos far from spectrum max_position (samples)" ) % j).str();
		    if ( ! m_bypassDebug )
		      Log.LogDebug ( "Peak_%d : gaussPos far from spectrum max_position (samples)", j );
                    result->PeakListDetectionStatus.push_back(status);
		  }
	      } // check gaussPos vs position of the max.
	    else
	      { //tolerance in Angstrom
                Float64 tolAngtsrom = 6;
                Float64 diffPos = fabs( gaussPos-spc.GetSpectralAxis()[max_index] );
                if( diffPos > tolAngtsrom )
		  {
                    toAdd = false;
                    std::string status = (boost::format( "Peak_%1% : gaussAmp far from spectrum max_value (Angstrom)" ) % j).str();
		    if ( ! m_bypassDebug )
		      Log.LogDebug ( "Peak_%d : gaussAmp far from spectrum max_value (Angstrom)", j );
                    result->PeakListDetectionStatus.push_back(status);
		  }
	      }
	  } // Check if gaussian fit is very different from peak itself

        // check type weak or strong
        Int32 force = 1; //weak by default
        Float64 ratioAmp = -1.0; //cut = -1.0 by default
        if( toAdd )
	  {
            ratioAmp = ComputeFluxes( spc, m_winsize, resPeaks[j] );
            if( ratioAmp<m_cut )
	      {
                toAdd = false;
                std::string status = (boost::format( "Peak_%d : ratioAmp<m_cut (%f<%f)" ) % j % ratioAmp % m_cut).str();
		if ( ! m_bypassDebug )
		  Log.LogDebug ( "Peak_%d : ratioAmp<m_cut (%f<%f)", j, ratioAmp, m_cut );
                result->PeakListDetectionStatus.push_back( status );
                // add this peak range to retest list
                retestPeaks.push_back( resPeaks[j] );
                retestGaussParams.push_back( SGaussParams( gaussPos, gaussAmp, gaussWidth ) );
	      }
	    else 
	      {
		if( ratioAmp>m_cut*m_strongcut )
		  {
		    force = 2; //strong
		  }
	      }
	  }

        if( toAdd )
	  {
            std::string status = (boost::format( "Peak_%1% : line detected successfully" ) % j).str();
            result->PeakListDetectionStatus.push_back( status );
            char buffer [64];
            sprintf( buffer,"detected_peak_%d", j );
            std::string peakName = buffer;
            result->RayCatalog.Add( CRay( peakName, gaussPos, m_type, std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM()), force , gaussAmp, gaussWidth, ratioAmp, gaussPosErr, gaussWidthErr, gaussAmpErr) );
	  }
    }

    // retest
    bool retest_flag = false;
    if( retestPeaks.size()>0 )
      {
        //CRayCatalog::TRayVector
        retest_flag = true;
      }
    while( retest_flag )
      {
        retest_flag = Retest( spectrum, *result, retestPeaks, retestGaussParams, result->RayCatalog.GetFilteredList( m_type, CRay::nForce_Strong ), m_winsize, m_cut );
      }

    return result;
}

/**
* Find the range that limits the gaussian fit according to input parameters.
* This function uses a maximum gaussian width param. in order to limit the size of the interval used for the fit.
*/
TInt32Range CLineDetection::LimitGaussianFitStartAndStop( Int32 i, const TInt32RangeList& peaksBorders, Int32 len, const CSpectrumSpectralAxis spectralAxis )
{
    Int32 fitStart = peaksBorders[i].GetBegin();
    Int32 fitStop = peaksBorders[i].GetEnd()+1;

    Float64 width = fitStop - fitStart;
    UInt32 center = fitStart + width/2;
    UInt32 start = std::max( 0, spectralAxis.GetIndexAtWaveLength( spectralAxis[center]-m_maxsize/2.0 ) );
    UInt32 stop = std::min( (Int32) len, spectralAxis.GetIndexAtWaveLength( spectralAxis[center]+m_maxsize/2.0 ) );
    Int32 maxwinsizeIndexes = stop-start;

    if( width>maxwinsizeIndexes )
      {
        fitStart = max( 0, (int) start );
        fitStop = min( len, (int) stop );
      }

    if( i>0 )
      {
        if( peaksBorders[i-1].GetEnd() > -1 )
	  {
	    fitStart = max( peaksBorders[i-1].GetEnd(), fitStart );
	  }
      }

    if( i<peaksBorders.size()-1 )
    {
        if( peaksBorders[i+1].GetBegin() > -1 )
            fitStop = min( peaksBorders[i+1].GetBegin(), fitStop );
    }

    return TInt32Range( fitStart, fitStop );
}

/**
 * Returns the ratio between max amplitude without continuum and noise.
 */
Float64 CLineDetection::ComputeFluxes( const CSpectrum& spectrum, Float64 winsize, TInt32Range range, TFloat64List mask, Float64 *maxFluxnoContinuum, Float64 *noise )
{
    const CSpectrum& spc = spectrum;
    const CSpectrumFluxAxis fluxAxis = spc.GetFluxAxis();
    const CSpectrumSpectralAxis specAxis = spc.GetSpectralAxis();

    //if no mask, then set it to 1
    if( mask.size()==0 )
      {
        mask.resize( spc.GetSampleCount() );
        for( Int32 k=0; k<spc.GetSampleCount(); k++ )
	  {
            mask[k] = 1;
	  }
      }

    Float64 maxValue = -DBL_MAX;
    Int32 maxIndex = -1;
    for( Int32 k=range.GetBegin(); k<=range.GetEnd(); k++ )
    {
        if( maxValue<fluxAxis[k] && mask[k]!=0 )
	  {
            maxValue = fluxAxis[k];
            maxIndex = k;
	  }
    }

    // strong/weak test to do
    const Float64* fluxData = fluxAxis.GetSamples();
    CMedian<Float64> medianProcessor;
    //reg. sampling, TAG: IRREGULAR_SAMPLING
    //Int32 windowSampleCount = winsize / spc.GetResolution();
    //int left = max(0, (Int32)(maxIndex-windowSampleCount/2.0+0.5) ) ;
    //int right = min((Int32)fluxAxis.GetSamplesCount()-1, (Int32)(maxIndex + windowSampleCount/2.0) )+1;
    //irreg. sampling
    int left = max( 0, specAxis.GetIndexAtWaveLength( specAxis[maxIndex]-winsize/2.0 ) );
    int right = min( (Int32)fluxAxis.GetSamplesCount(), specAxis.GetIndexAtWaveLength( specAxis[maxIndex]+winsize/2.0 )+1 );

    left = max(range.GetBegin(), left);
    right = min(range.GetEnd()+1, right);

    Float64 *fluxMasked = (Float64*)malloc( fluxAxis.GetSamplesCount()*sizeof( Float64 ) );
    int n=0;
    for( int i=left; i<right; i++ )
      {
        if( mask[i]!=0 )
	  {
            fluxMasked[n] = fluxData[i];
            n++;
	  }
      }
    Float64 med = medianProcessor.Find( fluxMasked, n );
    Float64 xmad = XMadFind( fluxMasked, n, med );
    free(fluxMasked);

    Float64 noise_win= xmad;
    // use noise spectrum
    const TFloat64List& error = fluxAxis.GetError().GetSamplesVector();
    //below operations can be moved to CSpectrumNoiseAxis
    if( !error.empty() )
      {
        // check if noise file has been loaded
        bool isNoiseOnes = true;
        for ( Int32 i=left; i<right; i++)
        {
            if( error[i]!=1.0 )
	          {
                isNoiseOnes = false;
                break;
	          }
        }

        if( !isNoiseOnes )
	      {
            Float64 mean_noise = 0.0;
            Int32 n_mean_noise = 0;
            for ( Int32 i=left; i<right; i++)
	          {
                if( mask[i]!=0 )
		            {
                    mean_noise += error[i];
                    n_mean_noise ++;
		            }
	          }
            if( n_mean_noise>0 )
	          {
                mean_noise /= n_mean_noise;
	          }
            // choose between noise mean or xmad
            if( mean_noise>xmad )
	          {
                noise_win = mean_noise;
	          }
            //noise_win = mean_noise; //debug
	      }

      }
    Float64 maxValue_no_continuum = maxValue - med;
    Float64 ratioAmp = maxValue_no_continuum/noise_win;

    // assign output pointers
    if( maxFluxnoContinuum != NULL )
      {
        *maxFluxnoContinuum = maxValue_no_continuum;
      }
    if( noise != NULL )
      {
        *noise = xmad;
      }

    return ratioAmp;
}

/**
 * Verifies the peaks are within the range of a strong peak.
 */
bool CLineDetection::Retest( const CSpectrum& spectrum, CLineDetectionResult& result, TInt32RangeList retestPeaks,  TGaussParamsList retestGaussParams, CRayCatalog::TRayVector strongLines, Int32 winsize, Float64 cut )
{
  if ( ! m_bypassDebug )
    Log.LogDebug( "Retest %d peaks, winsize = %d, strongLines.size() = %d.", retestPeaks.size(), winsize, strongLines.size() );
  DebugAssert( retestPeaks.size() == retestPeaks.size() );

  TGaussParamsList selectedgaussparams;
  TInt32RangeList selectedretestPeaks;

  const CSpectrumSpectralAxis wavesAxis = spectrum.GetSpectralAxis();
  // check if the retest peaks center are in the range of a strong peak
  for( int k=0; k<retestPeaks.size(); k++ )
      {
        Float64 start = wavesAxis[retestPeaks[k].GetBegin()];
        Float64 stop = wavesAxis[retestPeaks[k].GetEnd()];
        Float64 center = (stop + start)/2.0;

        for( int l=0; l<strongLines.size(); l++ )
	  {
            if( strongLines[l].GetPosition()-winsize/2.0<center && strongLines[l].GetPosition()+winsize/2.0 > center )
	      {
		if ( ! m_bypassDebug )
		  Log.LogDebug ( "The strongLine[%d].GetPosition() == %f is within winsize centered on %f.", l, strongLines[l].GetPosition(), center );
                selectedretestPeaks.push_back( retestPeaks[k] );
                selectedgaussparams.push_back( retestGaussParams[k] );
	      }
	    else
	      {
		if ( ! m_bypassDebug )
		  Log.LogDebug ( "The strongLine[%d].GetPosition() == %f is not within winsize centered on %f.", l, strongLines[l].GetPosition(), center );
	      }
	  }
      }

    if( selectedretestPeaks.size()<1 )
      {
	if ( ! m_bypassDebug )
	  Log.LogDebug ( "No retestPeaks were selected." );
        return false;
      }

    // remove selected retestPeaks: EZELFind ln 1102
    // not implemented, did not seem useful

    bool continue_retest;
    continue_retest = RemoveStrongFromSpectra( spectrum, result, strongLines, selectedretestPeaks, selectedgaussparams, winsize, cut );

    return continue_retest;
}

/**
 * Finds the intervals containing strong lines and removes the contents of those intervals from the spectrum.
 */
bool CLineDetection::RemoveStrongFromSpectra(const CSpectrum& spectrum, CLineDetectionResult& result,  CRayCatalog::TRayVector strongLines, TInt32RangeList selectedretestPeaks, TGaussParamsList selectedgaussparams, Float64 winsize, Float64 cut)
{

    const CSpectrumSpectralAxis wavesAxis = spectrum.GetSpectralAxis();
    // check intersection between stronglines and peaks selected
    TInt32RangeList toExclude;
    for( int l=0; l<strongLines.size(); l++ )
      {
        Float64 inf = strongLines[l].GetPosition() - FWHM_FACTOR*strongLines[l].GetCut();
        Float64 sup = strongLines[l].GetPosition() + FWHM_FACTOR*strongLines[l].GetCut();
        toExclude.push_back( TInt32Range( inf, sup ) );
    }
    // build toExclude to be non-intersecting stronglines
    for( int k=0; k<selectedretestPeaks.size(); k++ )
      {
        for( int l=0; l<toExclude.size(); l++ )
	  {
            TInt32Range line = TInt32Range( wavesAxis[selectedretestPeaks[k].GetBegin()], wavesAxis[selectedretestPeaks[k].GetEnd()] );
            TInt32Range strong = toExclude[l];
            if( line.GetEnd()>strong.GetBegin() && line.GetBegin()<strong.GetEnd() )
	      {
                if( strong.GetBegin() > line.GetBegin() )
		  {
                    toExclude[l].SetBegin( line.GetEnd() );
		  }
		else
		  {
                    toExclude[l].SetEnd( line.GetBegin() );
		  }
	      }
	  }
      }

    //create the mask on the strong peaks
    TFloat64List mask;
    mask.resize(spectrum.GetSampleCount());
    for( Int32 k=0; k<spectrum.GetSampleCount(); k++ )
    {
        mask[k] = 1;
    }
    for( Int32 k=0; k<wavesAxis.GetSamplesCount(); k++ )
    {
        for( int l=0; l<toExclude.size(); l++ )
	  {
            if( toExclude[l].GetBegin()<=wavesAxis[k] && toExclude[l].GetEnd()>=wavesAxis[k] )
	      {
                mask[k] = 0;
                break;
	      }
	  }
    }

    // create reduced spectra
    const CSpectrum reducedSpectrum( spectrum, mask );
    TFloat64List reducedindexesMap;
    for( Int32 k=0; k<spectrum.GetSampleCount(); k++ )
    {
        if( k==0 )
	  {
            reducedindexesMap.push_back( 0 );
	  }
	else
	  {
            if( mask[k] != 0 )
	      {
                reducedindexesMap.push_back( reducedindexesMap[k-1] + 1 );
	      }
	    else
	      {
                reducedindexesMap.push_back( reducedindexesMap[k-1] );
	      }

	  }
    }

    for( int k=0; k<selectedretestPeaks.size(); k++ )
      {
        TInt32Range reducedrange( (int)reducedindexesMap[selectedretestPeaks[k].GetBegin()], (int)reducedindexesMap[selectedretestPeaks[k].GetEnd()] );
        //Float64 ratioAmp = ComputeFluxes(spectrum, winsize, selectedretestPeaks[k], mask);
        Float64 ratioAmp = ComputeFluxes( reducedSpectrum, winsize, reducedrange );
        if( ratioAmp > cut )
	  {
            Float64 force = CRay::nForce_Weak;
            char buffer [64];
            sprintf( buffer, "detected_retested_peak_%d", k );
            std::string peakName = buffer;

            result.RayCatalog.Add( CRay( peakName, selectedgaussparams[k].Pos, m_type, std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM()), force , selectedgaussparams[k].Amp, selectedgaussparams[k].Width, ratioAmp) );
            result.RayCatalog.Sort();
        }
    }

    return false;
}

/**
 * Returns the median value of the residue of x with regards to median.
 */
Float64 CLineDetection::XMadFind( const Float64* x, Int32 n, Float64 median )
{
    std::vector<Float64> xdata;
    Float64 xmadm = 0.0;

    xdata.resize( n );

    for( Int32 i=0; i<n; i++ )
    {
        xdata[i] = fabs( x[i]-median );
    }

    std::sort(xdata.begin(), xdata.end());

    if ( n & 1 )
    {
        // n is odd
        UInt32 i1 = n >> 1; // i.e. int(n/2)
        xmadm = xdata[i1];
    }
    else
    {
        // n is even
        UInt32 i1 = n/2 -1;
        UInt32 i2 = n/2;
        xmadm = 0.5*(xdata[i1]+xdata[i2]);
    }

    return xmadm;
}

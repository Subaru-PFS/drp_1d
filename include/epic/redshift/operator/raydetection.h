#ifndef _REDSHIFT_OPERATOR_RAYDETECTION_
#define _REDSHIFT_OPERATOR_RAYDETECTION_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/spectrum/spectralaxis.h>

namespace NSEpic
{

class CSpectrum;
class CLineDetectionResult;

/**
 * \ingroup Redshift
 */
class CLineDetection
{
 public:
  struct SGaussParams
  {
    SGaussParams( )
    {
      Pos = 0.0;
      Width = 0.0;
      Amp = 0.0;
    }
    SGaussParams( Float64 pos, Float64 width, Float64 amp )
    {
      Pos = pos;
      Width = width;
      Amp = amp;
    }
    Float64 Pos;
    Float64 Width;
    Float64 Amp;
  };

  typedef std::vector<SGaussParams> TGaussParamsList;

  CLineDetection( Int32 type=CRay::nType_Emission, Float64 cut=5.0, Float64 strongcut=2.0, Float64 winsize=250, Float64 minsize=3, Float64 maxsize=70, bool disableFitQualityCheck=false );
  virtual ~CLineDetection();

  std::shared_ptr<const CLineDetectionResult> Compute( const CSpectrum& spectrum, const TLambdaRange& lambdaRange, const TInt32RangeList& resPeaks, const TInt32RangeList& resPeaksEnlarged );

  Float64 FWHM_FACTOR;
  
  Float64 ComputeFluxes( const CSpectrum& spectrum, Float64 winsize, TInt32Range range, TFloat64List mask=TFloat64List(),Float64* maxFluxnoContinuum=NULL, Float64* noise=NULL );
  
 private:
  Int32 m_type;
  
  Float64 m_winsize;
  Float64 m_minsize;
  Float64 m_maxsize;
  
  Float64 m_cut;
  Float64 m_strongcut;
  
  bool m_disableFitQualityCheck;
  
  TInt32Range LimitGaussianFitStartAndStop( Int32 i, const TInt32RangeList& peaksBorders, Int32 len, const CSpectrumSpectralAxis spectralAxis );
  
  bool Retest( const CSpectrum &spectrum, CLineDetectionResult& result, TInt32RangeList retestPeaks,  TGaussParamsList retestGaussParams, CRayCatalog::TRayVector strongLines, Int32 winsize, Float64 cut );
  bool RemoveStrongFromSpectra( const CSpectrum &spectrum, CLineDetectionResult& result, CRayCatalog::TRayVector strongLines, TInt32RangeList selectedretestPeaks, TGaussParamsList selectedgaussparams, Float64 winsize, Float64 cut );
  Float64 XMadFind( const Float64* x, Int32 n, Float64 median );

  // Log
  Bool m_bypassDebug; // If true, debug messages are ignored even if --verbose has been set.
};
}

#endif // _REDSHIFT_OPERATOR_RAYDETECTION_

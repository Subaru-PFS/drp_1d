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
#ifndef _REDSHIFT_OPERATOR_RAYDETECTION_
#define _REDSHIFT_OPERATOR_RAYDETECTION_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/range.h"

#include "RedshiftLibrary/operator/raydetectionresult.h"
#include "RedshiftLibrary/spectrum/spectralaxis.h"
//#include "RedshiftLibrary/ray/lineprofile.h"

#include "RedshiftLibrary/ray/lineprofileASYM.h"

namespace NSEpic
{

class CSpectrum;
class CLineDetectionResult;
class CLineProfileSYM;


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

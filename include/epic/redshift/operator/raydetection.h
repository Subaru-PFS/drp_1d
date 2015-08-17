#ifndef _REDSHIFT_OPERATOR_RAYDETECTION_
#define _REDSHIFT_OPERATOR_RAYDETECTION_

#include <epic/core/common/managedobject.h>
#include <epic/core/common/datatypes.h>
#include <epic/core/common/range.h>

#include <epic/redshift/operator/raydetectionresult.h>

namespace NSEpic
{

class CSpectrum;
class CRayDetectionResult;

class CRayDetection : public CManagedObject
{
    DEFINE_MANAGED_OBJECT( CRayDetection )

    public:

    struct SGaussParams
    {
        SGaussParams( )
        {
            Pos = 0.0;
            Width = 0.0;
            Amp = 0.0;
        }

        SGaussParams( Float64 pos, Float64 width, Float64 amp)
        {
            Pos = pos;
            Width = width;
            Amp = amp;
        }
    Float64 Pos;
    Float64 Width;
    Float64 Amp;
    };

    typedef std::vector<SGaussParams>   TGaussParamsList;

    CRayDetection(Float64 cut=5.0, Float64 strongcut=2.0, Float64 winsize=250, Float64 minsize=3, Float64 maxsize=70);
    virtual ~CRayDetection();

    const CRayDetectionResult* Compute(const CSpectrum& spectrum, const TLambdaRange& lambdaRange, const TInt32RangeList& resPeaks, const TInt32RangeList& resPeaksEnlarged);


    Float64 FWHM_FACTOR;

private:

    Float64 m_winsize;
    Float64 m_minsize;
    Float64 m_maxsize;

    Float64 m_cut;
    Float64 m_strongcut;


    Float64 ComputeFluxes(const CSpectrum& spectrum, Float64 winsize, TInt32Range range, TFloat64List mask=TFloat64List());
    bool Retest( const CSpectrum &spectrum, CRayDetectionResult* result, TInt32RangeList retestPeaks,  TGaussParamsList retestGaussParams, CRayCatalog::TRayVector strongLines, Int32 winsize, Float64 cut);
    bool RemoveStrongFromSpectra(const CSpectrum &spectrum, CRayDetectionResult* result, CRayCatalog::TRayVector strongLines, TInt32RangeList selectedretestPeaks, TGaussParamsList selectedgaussparams, Float64 winsize, Float64 cut);
    Float64 XMadFind( const Float64* x, Int32 n, Float64 median );

};

}

#endif

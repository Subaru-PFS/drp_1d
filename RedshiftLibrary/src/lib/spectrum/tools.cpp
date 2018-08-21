#include <RedshiftLibrary/spectrum/tools.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/spectrum/axis.h>
#include <RedshiftLibrary/common/mask.h>

#include <math.h>

using namespace NSEpic;

CSpectrumTools::CSpectrumTools()
{
}

CSpectrumTools::~CSpectrumTools()
{
}


void CSpectrumTools::Interpolate( const CSpectrumAxis& axisXorg, const CSpectrumAxis& axisYorg, Int32 offsetOrg, Int32 nOrg,
                                const CSpectrumAxis& axisXint, CSpectrumAxis& axisYint, CMask& mask )
{
    Int32 nInt = axisXint.GetSamplesCount();
    Int32 i; // move over original axis
    Int32 j; // move over interpolated axis

    const Float64* Xorg = axisXorg.GetSamples() + offsetOrg;
    const Float64* Yorg = axisYorg.GetSamples() + offsetOrg;
    const Float64* Xint = axisXint.GetSamples();
    Float64* Yint = axisYint.GetSamples();

    // While original interpolated axis doesn't overlap original axis,
    // mask this area
    for( i=0,j=0; j < nInt && Xint[j] < Xorg[i]; j++ )
    {
        mask[j]=0;
    }

    // For each sample in original axis
    for( i=0; i<nOrg-1; i++ )
    {
        // While spectrum points are between 2 template points
        while( j<nInt && Xint[j]<Xorg[i+1] )
        {
            // Perform linear interpolation of template flux value
            Float64 t = ( Xint[j]-Xorg[i] ) / ( Xorg[i+1]-Xorg[i] );

            TFloat64Range r = TFloat64Range( Yorg[i], Yorg[i+1] );
            Yint[j]= r.Interpolate( t );
            mask[j]=1;
            j++;
        }
    }

    for( j ;j<nInt; j++ )
    {
        //DebugAssert( Xint[j]>Xorg[Norg-1] );

        if (Xint[j]>Xorg[nOrg-1])
        {
            mask[j]=0;
        }
        else
        {
            mask[j]=1;
            Yint[j]=Yorg[j]; // No interpolation performed
        }

    }
}

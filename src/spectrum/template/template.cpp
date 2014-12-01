#include <epic/redshift/spectrum/template/template.h>

#include <epic/redshift/common/mask.h>

using namespace NSEpic;
using namespace std;

CTemplate::CTemplate() :
    m_Category( nCategory_None )
{

}

CTemplate::CTemplate( const char* name, ECategory category ) :
    m_Category( category ),
    m_Name( name )
{

}

CTemplate::~CTemplate()
{

}

Void CTemplate::Interpolate( const CTemplate& tpl, const CSpectrumSpectralAxis& targetSpectralAxis )
{
    GetFluxAxis().SetSize( targetSpectralAxis.GetSamplesCount() );
    GetSpectralAxis().CopyFrom( targetSpectralAxis );

    Int32 nOrg = tpl.GetSpectralAxis().GetSamplesCount();
    Int32 nInt = targetSpectralAxis.GetSamplesCount();
    Int32 i; // move over original axis
    Int32 j; // move over interpolated axis

    const Float64* Xorg = tpl.GetSpectralAxis().GetSamples() ;
    const Float64* Yorg = tpl.GetFluxAxis().GetSamples() ;
    Float64* Xint = GetSpectralAxis().GetSamples();
    Float64* Yint = GetFluxAxis().GetSamples();

    // While original interpolated axis doesn't overlap original axis,
    // mask this area
    for( i=0,j=0; j < nInt && Xint[j] < Xorg[i]; j++ )
    {
        Yint[j] = -1;
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
            j++;
        }
    }

    for( j ;j<nInt; j++ )
    {
        //DebugAssert( Xint[j]>Xorg[Norg-1] );

        if (Xint[j]>Xorg[nOrg-1])
        {
            Yint[j] = -1;
        }
        else
        {
            Yint[j]=Yorg[j]; // No interpolation performed
        }

    }
}

const char* CTemplate::GetCategoryName( ECategory cat )
{
    static const char* emission = "emission";
    static const char* galaxy = "galaxy";
    static const char* star = "star";
    static const char* qso = "qso";
    static const char* none = "none";

    if( cat == nCategory_Emission )
        return emission;
    else if( cat == nCategory_Galaxy )
        return galaxy;
    else if( cat == nCategory_Star )
        return star;
    else if( cat == nCategory_Qso )
        return qso;

    return none;
}

const std::string& CTemplate::GetName() const
{
    return m_Name;
}

CTemplate::ECategory CTemplate::GetCategory() const
{
    return m_Category;
}


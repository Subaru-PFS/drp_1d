#include <epic/redshift/method/linemodelsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/linemodelresult.h>
#include <stdio.h>
#include <float.h>
#include <math.h>


using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CLineModelSolveResult )

CLineModelSolveResult::CLineModelSolveResult()
{

}

CLineModelSolveResult::~CLineModelSolveResult()
{

}

Void CLineModelSolveResult::Save( const COperatorResultStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshiftBayes( store, redshift, merit );

    stream <<  "#Redshifts\tMerit\tTemplate"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << std::endl;


}


Void CLineModelSolveResult::SaveLine( const COperatorResultStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;

    GetBestRedshiftBayes( store, redshift, merit );

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "LineModelSolve" << std::endl;
}

Bool CLineModelSolveResult::GetBestRedshift( const COperatorResultStore& store, Float64& redshift, Float64& merit ) const
{

    std::string scope = store.GetScope( this ) + "linemodelsolve.linemodel";
    const CLineModelResult* results = (CLineModelResult*)store.GetGlobalResult(scope.c_str());


    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;

    if(results){
        for( Int32 i=0; i<results->ChiSquare.size(); i++ )
        {
            if( results->ChiSquare[i] < tmpMerit )
            {
                tmpMerit = results->ChiSquare[i];
                tmpRedshift = results->Redshifts[i];
            }
        }

    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;

}

Bool CLineModelSolveResult::GetBestRedshiftBayes( const COperatorResultStore& store, Float64& redshift, Float64& merit ) const
{

    std::string scope = store.GetScope( this ) + "linemodelsolve.linemodel";
    CLineModelResult* results = (CLineModelResult*)store.GetGlobalResult(scope.c_str());

    Float64 tmpMerit = DBL_MIN ;
    Float64 tmpRedshift = 0.0;

//    if(results){
//        //prepare p
//        Float64 minp = DBL_MIN;
//        CSpectrum pspc;
//        CSpectrumAxis& spcFluxAxis = pspc.GetFluxAxis();
//        spcFluxAxis.SetSize( results->Redshifts.size() );
//        CSpectrumAxis& spcSpectralAxis = pspc.GetSpectralAxis();
//        spcSpectralAxis.SetSize( results->Redshifts.size()  );
//        for( Int32 i2=0; i2<results->Redshifts.size(); i2++ )
//        {
//            spcFluxAxis[i2] = exp(-results->ChiSquare[i2]/2.0);
//            spcSpectralAxis[i2] = results->Redshifts[i2];
//            if( minp > spcFluxAxis[i2]){
//                minp = spcFluxAxis[i2];
//            }
//        }

//        Float64 winsize = 0.002;
//        for( Int32 i=0; i<results->Extrema.size(); i++ )
//        {
//            //find iz, izmin and izmax
//            Int32 izmin= -1;
//            Int32 iz= -1;
//            Int32 izmax= -1;
//            for( Int32 i2=0; i2<results->Redshifts.size(); i2++ )
//            {
//                if(iz == -1 && (results->Extrema[i]) <= results->Redshifts[i2]){
//                    iz = i2;
//                }
//                if(izmin == -1 && (results->Extrema[i] - winsize/2.0) <= results->Redshifts[i2]){
//                    izmin = i2;
//                }
//                if(izmax == -1 && (results->Extrema[i] + winsize/2.0) <= results->Redshifts[i2]){
//                    izmax = i2;
//                    break;
//                }
//            }

//            CGaussianFitSimple fitter;
//            CGaussianFitSimple::EStatus status = fitter.Compute( pspc, TInt32Range( izmin, izmax ) );
//            if(status!=NSEpic::CGaussianFitSimple::nStatus_Success){
//                continue;
//            }

//            Float64 gaussAmp;
//            Float64 gaussPos;
//            Float64 gaussWidth;
//            fitter.GetResults( gaussAmp, gaussPos, gaussWidth );
//            Float64 gaussAmpErr;
//            Float64 gaussPosErr;
//            Float64 gaussWidthErr;
//            fitter.GetResultsError( gaussAmpErr, gaussPosErr, gaussWidthErr );

//            Float64 area = gaussAmp*(gaussWidth*sqrt(2.0))*sqrt(2.0*3.141592654);
//            results->Area[i] = area;
//            if( area > tmpMerit ){
//                tmpMerit = area;
//                tmpRedshift = gaussPos;
//            }

//        }
//    }

    if(results){
        for( Int32 i=0; i<results->Area.size(); i++ )
        {
            if( results->Area[i] > tmpMerit )
            {
                tmpMerit = results->Area[i];
                tmpRedshift = results->Extrema[i];
            }
        }

    }

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
}




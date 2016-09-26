#include <epic/redshift/method/dtreecsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>

#include <epic/redshift/operator/linemodelresult.h>

#include <epic/core/log/log.h>

#include <float.h>

using namespace NSEpic;


CDTreeCSolveResult::CDTreeCSolveResult()
{

}

CDTreeCSolveResult::~CDTreeCSolveResult()
{

}

Void CDTreeCSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    std::string dtreepath;

    GetBestRedshift( store, redshift, merit, dtreepath );

    stream <<  "#Redshifts\tMerit\tTemplate"<< std::endl;

    stream  << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                   << "dtreec_" << dtreepath.c_str() << std::endl;

}

Void CDTreeCSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    std::string dtreepath;

    GetBestRedshift( store, redshift, merit, dtreepath );

    stream  << store.GetSpectrumName() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "dtreec_" << dtreepath.c_str() << std::endl;



    Log.LogInfo( "DecisionalTreeC Solution: best z found = %.5f", redshift);
}


Bool CDTreeCSolveResult::GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit , std::string &dtreepath) const
{
    //*
    // combined merit curve minimum
    std::string scope_lincomb = store.GetScope( *this ) + "dtreeCsolve.resultdtreeCCombined";

    auto _results = store.GetGlobalResult(scope_lincomb.c_str());
    auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>(_results.lock());



    //estimate the combined chisquare result
    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = -1.0;
    for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
    {
        Float64 post = meritResult->ChiSquare[i];
        if( post < tmpMerit )
        {
            tmpMerit = post;
            tmpRedshift = meritResult->Redshifts[i];
        }
    }


    redshift = tmpRedshift;
    merit = tmpMerit;
    dtreepath = "0.0";
    return true;
    //*/
}

///**
// * \brief Searches all the extrema results for the lowest, using the StrongELSnrPrior to correct the initial combined merit curve
// **/
//Bool CDTreeCSolveResult::GetBestRedshiftWithStrongELSnrPrior( const CDataStore& store, Float64& redshift, Float64& merit, std::string &dtreepath ) const
//{
//    // combined merit curve minimum
//    std::string scope_lincomb = store.GetScope( *this ) + "dtreeCsolve.resultdtreeCCombined";
//    auto _results_lincomb = store.GetGlobalResult(scope_lincomb.c_str());
//    auto meritResult_lincomb = std::dynamic_pointer_cast<const CChisquareResult>(_results_lincomb.lock());

//    std::string scope = store.GetScope( *this ) + "dtreeCsolve.linemodel";
//    auto results = store.GetGlobalResult( scope.c_str() );

//    Float64 tmpMerit = DBL_MAX ;
//    Float64 tmpRedshift = 0.0;

//    if(!results.expired())
//    {
//        auto lineModelResult = std::dynamic_pointer_cast<const CLineModelResult>( results.lock() );
//        for( Int32 i=0; i<lineModelResult->Extrema.size(); i++ )
//        {
//            //find lincomb corresponding index
//            Float64 idxLinComb = -1.0;
//            Float64 meritLinComb = DBL_MAX;
//            for( Int32 ilc=0; ilc<meritResult_lincomb->ChiSquare.size(); ilc++ )
//            {
//                if( lineModelResult->Extrema[i] == meritResult_lincomb->Redshifts[ilc])
//                {
//                    idxLinComb=ilc;
//                    meritLinComb = meritResult_lincomb->ChiSquare[ilc];
//                }
//            }

//            Float64 coeff =10.0;
//            Float64 correctedMerit = meritLinComb-coeff*lineModelResult->StrongELSNR[i];
//            if( correctedMerit < tmpMerit )
//            {
//                tmpMerit = correctedMerit;
//                tmpRedshift = lineModelResult->Extrema[i];
//            }
//        }
//    }

//    redshift = tmpRedshift;
//    merit = tmpMerit;
//    dtreepath = "0.0";
//    return true;
//}


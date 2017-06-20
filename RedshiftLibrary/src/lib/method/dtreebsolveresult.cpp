#include <RedshiftLibrary/method/dtreebsolveresult.h>

#include <RedshiftLibrary/processflow/context.h>
#include <RedshiftLibrary/operator/chisquareresult.h>
#include <RedshiftLibrary/operator/correlationresult.h>

#include <RedshiftLibrary/operator/linemodelresult.h>

#include <RedshiftLibrary/log/log.h>

#include <float.h>

using namespace NSEpic;


CDTreeBSolveResult::CDTreeBSolveResult()
{

}

CDTreeBSolveResult::~CDTreeBSolveResult()
{

}

Void CDTreeBSolveResult::Save( const CDataStore& store, std::ostream& stream ) const
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
                   << "dtreeb_" << dtreepath.c_str() << std::endl;

}

Void CDTreeBSolveResult::SaveLine( const CDataStore& store, std::ostream& stream ) const
{
    Float64 redshift;
    Float64 merit;
    std::string tplName;
    std::string dtreepath;

    GetBestRedshift( store, redshift, merit, dtreepath );

    stream  << store.GetSpectrumName() << "\t"
            << store.GetProcessingID() << "\t"
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "dtreeb_" << dtreepath.c_str() << std::endl;



    Log.LogInfo( "DecisionalTreeB Solution: best z found = %.5f", redshift);
}


Bool CDTreeBSolveResult::GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit , std::string &dtreepath) const
{
    //***********************************************************
    //first test, NStrong lines>3, keep linemodel
    /*
    std::string scope = store.GetScope( *this ) + "dtreeBsolve.linemodel";
    auto results = std::dynamic_pointer_cast<const CLineModelResult>( store.GetGlobalResult(scope.c_str()).lock() );

    Float64 SNRcutThresStrongLines = 5.0;
    Float64 FITcutThresStrongLines = 5.0;

    Int32 NStrongLinesThres = 3;
    std::vector<Int32> moreThanNStrongLinesSols;
    if(results){
        Int32 maxNStrongId=0;
        Int32 maxNStrong = 0;
        for( Int32 iE=0; iE<results->Extrema.size(); iE++ )
        {
//            if(!results->IsLocalExtrema[iE]){
//                continue;
//            }


            //print for debug
            //Int32 nStrongTmp = results->GetNLinesOverCutThreshold(iE, 5.0, 1.0);
            //Log.LogInfo( "dtreeBsolve nvalid strong test: z= %.4f\tnstrong= %d", results->Extrema[iE], nStrongTmp);


            Int32 nStrong = results->GetNLinesOverCutThreshold(iE, SNRcutThresStrongLines, FITcutThresStrongLines);
            if(nStrong>=NStrongLinesThres){
                moreThanNStrongLinesSols.push_back(iE);
                Log.LogInfo( "dtreeBsolve gbr: z= %.4f\tnstrong= %d", results->Extrema[iE], nStrong);
                if(maxNStrong<nStrong){
                    maxNStrongId = iE;
                    maxNStrong = nStrong;
                }
            }

        }
        Log.LogInfo( "\n");
        if(moreThanNStrongLinesSols.size()>=1){
            redshift = results->Extrema[maxNStrongId];
            merit = 1;
            Log.LogInfo( "dtreeBsolve : Found N valid Strong lines result, z=%f", redshift);
            dtreepath = "1.1";
            return true;
        }else{
            redshift = 0;
            merit = 0;
        }
    }
    //*/


    //*
    //***********************************************************
    // Next Test: linear combination chisquare
    std::string scope_lincomb = store.GetScope( *this ) + "dtreeBsolve.resultdtreeBCombined";

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
    dtreepath = "3.1";
    return true;
    //*/

}


#include <epic/redshift/method/dtreebsolveresult.h>

#include <epic/redshift/processflow/context.h>
#include <epic/redshift/operator/chisquareresult.h>
#include <epic/redshift/operator/correlationresult.h>

#include <epic/redshift/operator/linemodelresult.h>

#include <epic/core/log/log.h>

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
                << redshift << "\t"
                << merit << "\t"
                << tplName << "\t"
                << "dtreeb_" << dtreepath.c_str() << std::endl;



    Log.LogInfo( "dtreeBsolve Solution: best z found = %.5f", redshift);
}


Bool CDTreeBSolveResult::GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit , std::string &dtreepath) const
{
    std::string scope = store.GetScope( *this ) + "dtreeBsolve.linemodel";
    auto results = std::dynamic_pointer_cast<const CLineModelResult>( store.GetGlobalResult(scope.c_str()).lock() );

    Float64 SNRcutThresStrongLines = 5.0;
    Float64 FITcutThresStrongLines = 5.0;

    //***********************************************************
    //first test, NStrong lines>3, keep linemodel
    //*
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

            /*
            //print for debug
            Int32 nStrongTmp = results->GetNLinesOverCutThreshold(iE, 5.0, 1.0);
            Log.LogInfo( "dtreeBsolve nvalid strong test: z= %.4f\tnstrong= %d", results->Extrema[iE], nStrongTmp);
            */

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


    /*
    //***********************************************************
    //second test, chi2 diff > 300 then keep the first extrema
    Float64 chi2diffThresLevel1 = 100;
    Float64 chi2diffThresLevel2 = 50;
    Float64 bestExtremaMerit = DBL_MAX;
    Float64 thres = 0.002;
    Int32 idxNextValid = 1;
    Int32 idxFirstValid = 1;
    Float64 firstz = -1;
    Int32 nExtrema = results->Extrema.size();
    for( Int32 iE=0; iE<nExtrema; iE++ )
    {
        if(firstz == -1 || bestExtremaMerit > results->GetExtremaMerit(iE)){
            firstz = results->Extrema[iE];
            bestExtremaMerit = results->GetExtremaMerit(iE);
            idxFirstValid = iE;
        }
    }
    Log.LogInfo( "dtreeBsolve : bestExtremaMerit, %f", bestExtremaMerit);

    Float64 nextExtremaMerit = DBL_MAX;
    for( Int32 iE=0; iE<nExtrema; iE++ )
    {

        Float64 thisz = results->Extrema[iE];

        Float64 diffz = std::abs(thisz-firstz);
        if( diffz > thres && nextExtremaMerit>results->GetExtremaMerit(iE)){
            idxNextValid = iE;
            nextExtremaMerit = results->GetExtremaMerit(iE);
        }
    }
    Log.LogInfo( "dtreeBsolve : nextExtremaMerit, %f", nextExtremaMerit);
    Float64 chi2diff = -(bestExtremaMerit - nextExtremaMerit);
    Int32 nStrongFirstExtremum = results->GetNLinesOverCutThreshold(idxFirstValid, 5.0, 1.5);
    Int32 nStrongNextExtremum = results->GetNLinesOverCutThreshold(idxNextValid, 5.0, 1.5);
    if( chi2diff > chi2diffThresLevel1 ){
        redshift = firstz;
        merit = bestExtremaMerit;
        Log.LogInfo( "dtreeBsolve : Found chi2diff result, z=%f", redshift);
        dtreepath = "2.1";
        return true;
    }else if(chi2diff > chi2diffThresLevel2 && nStrongFirstExtremum>=2 && nStrongNextExtremum<nStrongFirstExtremum){
        redshift = firstz;
        merit = bestExtremaMerit;
        Log.LogInfo( "dtreeBsolve : Found chi2diff result, with nstrongextr=%d, z=%f", nStrongFirstExtremum, redshift);
        dtreepath = "2.2";
        return true;
    }
//    else{
//        redshift = 0;
//        merit = 0;
//        return true;
//    }
    //*/

    //*
    //***********************************************************
    //retrieve chi2nc values
    std::vector<Float64> w_chi2nc;
    TFloat64List chi2nc;
    TFloat64List znc;
    Float64 minchi2nc= DBL_MAX;

    chi2nc = GetBestRedshiftChi2List(store, "chisquare_nocontinuum", minchi2nc, znc);

    //Log.LogInfo( "dtreeBsolve : znc size=%d", znc.size());
    for( Int32 i=0; i<chi2nc.size(); i++ )
    {
        w_chi2nc.push_back(chi2nc[i]/minchi2nc);
    }

    /*
    //***********************************************************
    //retrieve chi2 raw values
    std::vector<Float64> w_chi2raw;
    TFloat64List chi2raw;
    TFloat64List zraw;
    Float64 minchi2raw= DBL_MAX;

    chi2raw = GetBestRedshiftChi2List(store, "chisquare", minchi2raw, zraw);
    for( Int32 i=0; i<chi2raw.size(); i++ )
    {
        w_chi2raw.push_back(chi2raw[i]/minchi2raw);
    }
    //*/


    //*
    //***********************************************************
    //retrieve chi2 continuum values
    std::vector<Float64> w_chi2continuum;
    TFloat64List chi2continuum;
    TFloat64List zcontinuum;
    Float64 minchi2continuum= DBL_MAX;

    chi2continuum = GetBestRedshiftChi2List(store, "chisquare_continuum", minchi2continuum, zcontinuum);
    for( Int32 i=0; i<chi2continuum.size(); i++ )
    {
        w_chi2continuum.push_back(chi2continuum[i]/minchi2continuum);
    }

    //*
    //***********************************************************
    //retrieve linemodel values
    TFloat64List chi2lm;

    for( Int32 i=0; i<znc.size(); i++ )
    {
        for( Int32 iall=0; iall<results->Redshifts.size(); iall++ )
        {
            if(results->Redshifts[iall] == znc[i]){
                chi2lm.push_back(results->ChiSquare[iall]);
                break;
            }
        }
    }
    //Log.LogInfo( "dtreeBsolve : chi2lm size=%d", chi2lm.size());

    /*
    // save all merits in a temp. txt file
    FILE* f = fopen( "dtreeb_merits_dbg.txt", "w+" );
    for( Int32 i=0; i<results->Redshifts.size(); i++ )
    {
        fprintf( f, "%i %f %f %f %f\n", i, results->Extrema[i], results->GetExtremaMerit(i)/((float)results->nSpcSamples), chi2nc[i], chi2continuum[i] );//*1e12);
    }
    fclose( f );
    //*/

    //*
    //***********************************************************
    // Next Test: chi2nc chi2cont linear combination
    //*
    Float64 lmCoeff = 75.0/((float)results->nSpcSamples);
    Float64 chi2ncCoeff = 100.0;
    Float64 chi2cCoeff = 50.0;
    //*/
    /*
    Float64 lmCoeff = 12.0/((float)results->nSpcSamples);
    Float64 chi2ncCoeff = 100.0;
    Float64 chi2cCoeff = 8.0;
    //*/

    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = -1.0;
    for( Int32 i=0; i<znc.size(); i++ )
    {
        Float64 post = lmCoeff*chi2lm[i] + chi2ncCoeff*chi2nc[i] + chi2cCoeff*chi2continuum[i];

        if( post < tmpMerit )
        {
            tmpMerit = post;
            tmpRedshift = znc[i];
        }
    }
    redshift = tmpRedshift;
    merit = tmpMerit;
    dtreepath = "3.1";
    return true;
    //*/

    /*
    //***********************************************************
    // Next Test: if chi2cont and linemodel agree
    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = -1.0;
//    if(w_chi2raw[0] == 1.0){
//        merit = results->GetExtremaMerit(idxNextValid);
//        redshift = results->Extrema[0];
//        return true;
//    }
    Float64 lmMeritRatioTol = 0.01;
    for( Int32 iE=0; iE<results->Extrema.size(); iE++ )
    {
        Float64 lmMeritRatio = results->GetExtremaMerit(iE)/results->GetExtremaMerit(0);
        if(w_chi2raw[iE] == 1.0 && lmMeritRatio<(1.0+lmMeritRatioTol)){
            tmpMerit = results->GetExtremaMerit(idxNextValid);
            tmpRedshift = results->Extrema[iE];
        }
    }
    if(tmpRedshift > -1.0){
        redshift = tmpRedshift;
        merit = tmpMerit;
        return true;
    }
    //*/


    /*
    //***********************************************************
    // Next Test: if chi2nc and linemodel agree
    tmpMerit = DBL_MAX ;
    tmpRedshift = -1.0;
//    if(w_chi2raw[0] == 1.0){
//        merit = results->GetExtremaMerit(idxNextValid);
//        redshift = results->Extrema[0];
//        return true;
//    }
    Float64 lmMeritRatioTol_nc = 0.01;
    for( Int32 iE=0; iE<results->Extrema.size(); iE++ )
    {
        Float64 lmMeritRatio = results->GetExtremaMerit(iE)/results->GetExtremaMerit(0);
        if(w_chi2nc[iE] == 1.0 && lmMeritRatio<(1.0+lmMeritRatioTol_nc)){
            tmpMerit = results->GetExtremaMerit(idxNextValid);
            tmpRedshift = results->Extrema[iE];
        }
    }
    if(tmpRedshift > -1.0){
        redshift = tmpRedshift;
        merit = tmpMerit;
        return true;
    }
    //*/

    /*
    //***********************************************************
    //then, use chi2 no continuum
    if(results){
        for( Int32 iE=0; iE<results->Extrema.size(); iE++ )
        {
//            Int32 i=0;
//            for ( UInt32 i2=0; i2<results->LineModelSolutions.size(); i2++)
//            {
//                if(results->Redshifts[i2] == results->Extrema[iE]){
//                    i = i2;
//                    break;
//                }
//            }


            std::string scopeStr;
            //get the no_continuum redshift result
            Float64 redshift_nc;
            Float64 merit_nc;
            std::string tplName_nc;
            //scopeStr = "chisquare";
            scopeStr = "chisquare_nocontinuum";
            //scopeStr = "chisquare_continuum";
            GetBestRedshiftChi2( store, scopeStr, results->Extrema[iE], redshift_nc, merit_nc, tplName_nc );

            //get the raw redshift result
            Float64 redshift_raw;
            Float64 merit_raw;
            std::string tplName_raw;
            scopeStr = "chisquare";
            GetBestRedshiftChi2( store, scopeStr, results->Extrema[iE], redshift_raw, merit_raw, tplName_raw );

            //get the continuum redshift result
            Float64 redshift_continuum;
            Float64 merit_continuum;
            std::string tplName_continuum;
            scopeStr = "chisquare_continuum";
            GetBestRedshiftChi2( store, scopeStr, results->Extrema[iE], redshift_continuum, merit_continuum, tplName_continuum );


            Int32 nStrong = results->GetNLinesOverCutThreshold(iE, cutThresStrongLines);
            Log.LogInfo( "dtreeBsolve gbr: z= %.4f\tnstrong= %d\tmerit_linemodel= %.5f\t merit_chi2nc= %.5f\tmerit_chi2= %.5f\tmerit_chi2cont= %.5f\t", results->Extrema[iE], nStrong, results->GetExtremaMerit(iE), merit_nc, merit_raw, merit_continuum);


            {
                //if(nStrong>=2){
                Float64 post = merit_nc;
                //Float64 post = results->GetExtremaMerit(iE);

                if( post < tmpMerit )
                {
                    tmpMerit = post;
                    tmpRedshift = results->Extrema[iE];
                }
            }


//            if(nStrong>=2){
//                Float64 post = results->Posterior[iE];

//                if( post < tmpMerit )
//                {
//                    tmpMerit = post;
//                    tmpRedshift = results->Extrema[iE];
//                }
//            }else{
//                redshift = 0;
//                merit = 0;
//                return true;
//            }
        }

    }
    Log.LogInfo( "\n");

    redshift = tmpRedshift;
    merit = tmpMerit;
    return true;
    //*/

}


Bool CDTreeBSolveResult::GetBestRedshiftChi2( const CDataStore& store, std::string scopeStr, Float64 targetz, Float64& redshift, Float64& merit, std::string& tplName ) const
{
    std::string scope = store.GetScope( *this ) + "dtreeBsolve.chisquare2solve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    Float64 tmpMerit = DBL_MAX ;
    Float64 tmpRedshift = 0.0;
    std::string tmpTplName;

    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        Float64 zthres = 1e-12;
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>((*it).second);

        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if(std::abs(targetz-meritResult->Redshifts[i])<zthres && tmpMerit>meritResult->ChiSquare[i]){
                tmpMerit = meritResult->ChiSquare[i];
                tmpRedshift = meritResult->Redshifts[i];
                tmpTplName = (*it).first;
            }
        }
    }


    if( tmpMerit < DBL_MAX )
    {
        redshift = tmpRedshift;
        merit = tmpMerit;
        tplName = tmpTplName;
        return true;
    }

    return false;

}

TFloat64List CDTreeBSolveResult::GetBestRedshiftChi2List( const CDataStore& store, std::string scopeStr,  Float64& minmerit, TFloat64List& zList) const
{
    std::string scope = store.GetScope( *this ) + "dtreeBsolve.chisquare2solve." + scopeStr.c_str();
    TOperatorResultMap meritResults = store.GetPerTemplateResult(scope.c_str());

    TFloat64List meritList;
    //init meritresults
    TOperatorResultMap::const_iterator it0 = meritResults.begin();
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult> ( (*it0).second );
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            meritList.push_back(DBL_MAX);
            zList.push_back(meritResult->Redshifts[i]);
        }
    }

    //find best merit for each tpl
    minmerit = DBL_MAX;
    for( TOperatorResultMap::const_iterator it = meritResults.begin(); it != meritResults.end(); it++ )
    {
        auto meritResult = std::dynamic_pointer_cast<const CChisquareResult>((*it).second);
        for( Int32 i=0; i<meritResult->ChiSquare.size(); i++ )
        {
            if( meritList[i] > meritResult->ChiSquare[i]){
                meritList[i] = meritResult->ChiSquare[i];
            }
            if( minmerit > meritResult->ChiSquare[i]){
                minmerit = meritResult->ChiSquare[i];
            }
        }
    }

    return meritList;

}

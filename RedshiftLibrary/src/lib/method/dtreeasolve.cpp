#include <RedshiftLibrary/method/dtreeasolve.h>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chicorr.h>
#include <RedshiftLibrary/method/blindsolveresult.h>

#include <RedshiftLibrary/method/blindsolve.h>
#include <RedshiftLibrary/method/blindsolveresult.h>
#include <RedshiftLibrary/operator/chisquare.h>

#include <RedshiftLibrary/operator/peakdetection.h>
#include <RedshiftLibrary/operator/peakdetectionresult.h>
#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raymatching.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>

#include <RedshiftLibrary/method/chisquare2solve.h>
#include <RedshiftLibrary/method/fullsolve.h>
#include <RedshiftLibrary/method/correlationsolve.h>
#include <RedshiftLibrary/method/linematchingsolve.h>
#include <RedshiftLibrary/method/dtreeasolve.h>
#include <RedshiftLibrary/method/dtreeasolveresult.h>

using namespace NSEpic;
using namespace std;


COperatorDTreeASolve::COperatorDTreeASolve(string calibrationPath)
{
    m_calibrationPath = calibrationPath;

    // Peak Detection
    m_winsize = 250.0;
    m_cut = 1.5;
    m_strongcut = 2.0;
    m_minsize = 3;
    m_maxsize = 70;
    m_enlargeRate = 2.0;

    // Line Matching
    m_minMatchNum = 1;
    m_tol = 0.002;


    //pfs TF overrides
    if(1)
    {
        // TF
        //m_winsize = 250.0;
        //m_cut = 150;
        //m_maxsize = 120;
        //m_enlargeRate = 1.0;
        //m_tol = 0.0025;

        // F + ErrF
        m_winsize = 250.0;
        m_cut = 1.5;
        m_maxsize = 120;
        m_enlargeRate = 2.0;
        m_tol = 0.0025;
    }
}

COperatorDTreeASolve::~COperatorDTreeASolve()
{

}

std::shared_ptr<const CDTreeASolveResult> COperatorDTreeASolve::Compute(CDataStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                                                        Int32 correlationExtremumCount, Float64 overlapThreshold )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( resultStore, "dtreeAsolve" );

    storeResult = Solve(resultStore, spc, spcWithoutCont,
                                       tplCatalog, tplCategoryList, restRayCatalog,
                                       lambdaRange, redshiftRange, redshiftStep,
                                       correlationExtremumCount, overlapThreshold );

    //storeResult = true;
    if( storeResult )
    {
        return std::shared_ptr<const CDTreeASolveResult>( new CDTreeASolveResult() );
    }

    return NULL;
}

Bool COperatorDTreeASolve::Solve(CDataStore &resultStore, const CSpectrum &spc, const CSpectrum &spcWithoutCont,
                                 const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList,
                                 const CRayCatalog &restRayCatalog, const TFloat64Range &lambdaRange,
                                 const TFloat64Range &redshiftRange, Float64 redshiftStep, Int32 correlationExtremumCount, Float64 overlapThreshold)
{

    TStringList   filteredTemplateCategoryList = getFilteredTplCategory( tplCategoryList, "emission" );

    CPeakDetection peakDetection(m_winsize, m_cut, 1, m_enlargeRate);
    auto peakDetectionResult = peakDetection.Compute( spc, lambdaRange);
    resultStore.StoreScopedGlobalResult( "peakdetection", peakDetectionResult );
    Log.LogInfo( "DTreeA - Peak Detection output: %d peaks found", peakDetectionResult->PeakList.size());

    // prepare the unused masks
    std::vector<CMask> maskList;

    // check Peak Detection results
    if(peakDetectionResult->PeakList.size()<1){
        Log.LogInfo( "DTreeA - No Peak found, switching to Chisquare Solve");
        resultStore.SetParam( "dtreepathnum", 1.1 );
        {
            //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquare2Solve chiSolve(m_calibrationPath);
            auto chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold, maskList );
            if( chisolveResult ) {
                resultStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
            }
            return true;
        }
    }

    // Line Detection
    CLineDetection lineDetection(CRay::nType_Emission, m_cut, m_strongcut, m_winsize, m_minsize, m_maxsize);
    auto lineDetectionResult = lineDetection.Compute( spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );
    resultStore.StoreScopedGlobalResult( "raycatalog", lineDetectionResult );
    Log.LogInfo( "DTreeA - Line Detection output: %d line(s) found", lineDetectionResult->RayCatalog.GetList().size());

    // check line Detection results
    Int32 nRaysDetected = lineDetectionResult->RayCatalog.GetList().size();
    if( nRaysDetected < 1){
        Log.LogInfo( "DTreeA - No line found, switching to Chisquare Solve");
        resultStore.SetParam( "dtreepathnum", 1.11 );
        {
            //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquare2Solve chiSolve(m_calibrationPath);
            auto chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold, maskList );
            if( chisolveResult ) {
                resultStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
            }
            return true;
        }
    }

    // Line Match
    Int32 typeFilter = CRay::nType_Emission;
    Int32 restForceFilter = -1;
    if(nRaysDetected == 1){
        restForceFilter =   CRay::nForce_Strong;
    }

    CRayMatching rayMatching;
    auto rayMatchingResult = rayMatching.Compute( lineDetectionResult->RayCatalog, restRayCatalog, redshiftRange, m_minMatchNum, m_tol, typeFilter, -1, restForceFilter);

    if(rayMatchingResult!=NULL){
        // Store matching results
        resultStore.StoreScopedGlobalResult( "raymatching", rayMatchingResult );

        //check line matching results
        if(rayMatchingResult->GetSolutionsListOverNumber(0).size()<1){
            Log.LogInfo( "DTreeA - Not match found [1], switching to Chisquare");
            resultStore.SetParam( "dtreepathnum", 3.1 );
            {  //chisquare
                // Create redshift initial list by spanning redshift acdross the given range, with the given delta
                TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
                DebugAssert( redshifts.size() > 0 );

                CMethodChisquare2Solve chiSolve(m_calibrationPath);
                auto chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                    tplCatalog, tplCategoryList,
                                                                                    lambdaRange, redshifts, overlapThreshold, maskList );
                if( chisolveResult ) {
                    resultStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
                }
                return true;
            }
        }
    }else{ //TBD: is that likely to happen ?
        Log.LogInfo( "DTreeA - Not match found  [0], switching to Chisquare");
        resultStore.SetParam( "dtreepathnum", 3.11 );
        {  //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquare2Solve chiSolve(m_calibrationPath);
            auto chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold, maskList );
            if( chisolveResult ) {
                resultStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
            }
            return true;
        }
    }

    Int32 matchNum = rayMatchingResult->GetMaxMatchingNumber();
    UInt32 nStrongPeaks = lineDetectionResult->RayCatalog.GetFilteredList(CRay::nType_Emission, CRay::nForce_Strong).size();

    //
    if(matchNum >= 1){
        resultStore.SetParam( "dtreepathnum", 4.1 );
        Log.LogInfo( "DTreeA - compute chisquare on redshift candidates from line matching" );
        rayMatchingResult->FilterWithRules(spc, lambdaRange, m_winsize);
        TFloat64List redshifts = rayMatchingResult->GetExtendedRedshiftCandidatesOverNumber(0, redshiftStep, 0.01);
        Log.LogInfo( "DTreeA - (n candidates = %d)", redshifts.size());
        { //chisolve with emission templates only
            CMethodChisquare2Solve chiSolve(m_calibrationPath);
            auto chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold, maskList );
            if( chisolveResult ) {
                resultStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
            }
            return true;
        }
    }

    Log.LogInfo( "DTreeA - no other path found than switching to Chisquare...");
    resultStore.SetParam( "dtreepathnum", 5.1 );
    {  //chisquare
        // Create redshift initial list by spanning redshift acdross the given range, with the given delta
        TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
        DebugAssert( redshifts.size() > 0 );

        CMethodChisquare2Solve chiSolve(m_calibrationPath);
        auto chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                            tplCatalog, tplCategoryList,
                                                                            lambdaRange, redshifts, overlapThreshold, maskList );
        if( chisolveResult ) {
            resultStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
        }
        return true;
    }
}


TStringList COperatorDTreeASolve::getFilteredTplCategory( const TStringList& tplCategoryListIn, const std::string& CategoryFilter)
{
    TStringList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<tplCategoryListIn.size(); i++ )
    {
        std::string category = tplCategoryListIn[i];
        if( category == "star" )
        {
        }
        else if(CategoryFilter == "all" || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    return filteredTemplateCategoryList;
}

#include <epic/redshift/method/dtreeasolve.h>

#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/operator/resultstore.h>

#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chicorr.h>
#include <epic/redshift/method/blindsolveresult.h>

#include <epic/redshift/method/blindsolve.h>
#include <epic/redshift/method/blindsolveresult.h>
#include <epic/redshift/operator/chisquare.h>

#include <epic/redshift/operator/peakdetection.h>
#include <epic/redshift/operator/peakdetectionresult.h>
#include <epic/redshift/operator/raydetection.h>
#include <epic/redshift/operator/raydetectionresult.h>
#include <epic/redshift/operator/raymatching.h>
#include <epic/redshift/operator/raymatchingresult.h>

#include <epic/redshift/method/chisquaresolve.h>
#include <epic/redshift/method/fullsolve.h>
#include <epic/redshift/method/correlationsolve.h>
#include <epic/redshift/method/linematchingsolve.h>
#include <epic/redshift/method/dtreeasolve.h>
#include <epic/redshift/method/dtreeasolveresult.h>

using namespace NSEpic;
using namespace std;

IMPLEMENT_MANAGED_OBJECT( COperatorDTreeASolve )

COperatorDTreeASolve::COperatorDTreeASolve()
{
    // Peak Detection
    m_winsize = 250.0;
    m_cut = 5.0;
    m_strongcut = 2.0;

    // Line Matching
    m_minMatchNum = 1;
    m_tol = 0.002;

    //dtree path
    m_dtreepathnum = -1.0;
}

COperatorDTreeASolve::~COperatorDTreeASolve()
{

}

const CDTreeASolveResult* COperatorDTreeASolve::Compute(COperatorResultStore& resultStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TTemplateCategoryList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                                                        Int32 correlationExtremumCount, Float64 overlapThreshold )
{
    Bool storeResult = false;

    COperatorResultStore::CAutoScope resultScope( resultStore, "dtreeAsolve" );

    storeResult = Solve(resultStore, spc, spcWithoutCont,
                                       tplCatalog, tplCategoryList, restRayCatalog,
                                       lambdaRange, redshiftRange, redshiftStep,
                                       correlationExtremumCount, overlapThreshold );

    //storeResult = true;
    if( storeResult )
    {
        CDTreeASolveResult*  SolveResult = new CDTreeASolveResult();
        return SolveResult;
    }

    return NULL;
}

Bool COperatorDTreeASolve::Solve(COperatorResultStore &resultStore, const CSpectrum &spc, const CSpectrum &spcWithoutCont, const CTemplateCatalog &tplCatalog, const TTemplateCategoryList &tplCategoryList, const CRayCatalog &restRayCatalog, const TFloat64Range &lambdaRange, const TFloat64Range &redshiftRange, Float64 redshiftStep, Int32 correlationExtremumCount, Float64 overlapThreshold)
{

    TTemplateCategoryList   filteredTemplateCategoryList = getFilteredTplCategory( tplCategoryList, CTemplate::nCategory_Emission);

    CPeakDetection peakDetection(m_winsize, m_cut);
    CConstRef<CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( spc, lambdaRange);
    resultStore.StoreGlobalResult( "peakdetection", *peakDetectionResult );
    Log.LogInfo( "DTreeA - Peak Detection output: %d peaks found", peakDetectionResult->PeakList.size());

    // check Peak Detection results
    if(peakDetectionResult->PeakList.size()<1){
        Log.LogInfo( "DTreeA - No Peak found, switching to Chisquare Solve");
        m_dtreepathnum = 1.1;
        {
            //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquareSolve chiSolve;
            CConstRef<CChisquareSolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold );
            if( chisolveResult ) {
                resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
            }
            return true;
        }
    }

    // Ray Detection
    CRayDetection rayDetection( m_cut, m_strongcut );
    CConstRef<CRayDetectionResult> rayDetectionResult = rayDetection.Compute( spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );
    resultStore.StoreGlobalResult( "raycatalog", *rayDetectionResult );
    Log.LogInfo( "DTreeA - Ray Detection output: %d ray(s) found", rayDetectionResult->RayCatalog.GetList().size());

    // check Ray Detection results
    Int32 nRaysDetected = rayDetectionResult->RayCatalog.GetList().size();
    if( nRaysDetected < 1){
        Log.LogInfo( "DTreeA - No ray found, switching to Chisquare Solve");
        m_dtreepathnum = 1.11;
        {
            //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquareSolve chiSolve;
            CConstRef<CChisquareSolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold );
            if( chisolveResult ) {
                resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
            }
            return true;
        }
    }

    // Ray Match
    CRayMatching rayMatching;
    CRef<CRayMatchingResult> rayMatchingResult = rayMatching.Compute(rayDetectionResult->RayCatalog, restRayCatalog, redshiftRange, m_minMatchNum, m_tol );
    if(rayMatchingResult!=NULL){
        // Store matching results
        resultStore.StoreGlobalResult( "raymatching", *rayMatchingResult );

        //check ray matching results
        if(rayMatchingResult->GetSolutionsListOverNumber(0).size()<1){
            Log.LogInfo( "DTreeA - Not match found [1], switching to Chisquare");
            m_dtreepathnum = 3.1;
            {  //chisquare
                // Create redshift initial list by spanning redshift acdross the given range, with the given delta
                TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
                DebugAssert( redshifts.size() > 0 );

                CMethodChisquareSolve chiSolve;
                CConstRef<CChisquareSolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                    tplCatalog, filteredTemplateCategoryList,
                                                                                    lambdaRange, redshifts, overlapThreshold );
                if( chisolveResult ) {
                    resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
                }
                return true;
            }
        }
    }else{
        Log.LogInfo( "DTreeA - Not match found  [0], switching to Chisquare");
        m_dtreepathnum = 3.11;
        {  //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquareSolve chiSolve;
            CConstRef<CChisquareSolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold );
            if( chisolveResult ) {
                resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
            }
            return true;
        }
    }

    Int32 matchNum = rayMatchingResult->GetMaxMatchingNumber();
    UInt32 nStrongPeaks = rayDetectionResult->RayCatalog.GetFilteredList(CRay::nType_Emission, CRay::nForce_Strong).size();

    // match num >= 2
    if(matchNum >= 2){
        m_dtreepathnum = 4.1;
        Log.LogInfo( "DTreeA - compute chisquare on redshift candidates from ray matching" );
        TFloat64List roundedRedshift = rayMatchingResult->GetRoundedRedshiftCandidatesOverNumber(2-1, redshiftStep);
        Log.LogInfo( "DTreeA - (n candidates = %d)", roundedRedshift.size());
        { //chisolve with emission templates only
            CMethodChisquareSolve chiSolve;
            CConstRef<CChisquareSolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, roundedRedshift, overlapThreshold );
            if( chisolveResult ) {
                resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
            }
            return true;
        }
    }

    Log.LogInfo( "DTreeA - no other path found than switching to Chisquare...");
    m_dtreepathnum = 3.12;
    {  //chisquare
        // Create redshift initial list by spanning redshift acdross the given range, with the given delta
        TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
        DebugAssert( redshifts.size() > 0 );

        CMethodChisquareSolve chiSolve;
        CConstRef<CChisquareSolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                            tplCatalog, filteredTemplateCategoryList,
                                                                            lambdaRange, redshifts, overlapThreshold );
        if( chisolveResult ) {
            resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
        }
        return true;
    }
}


TTemplateCategoryList COperatorDTreeASolve::getFilteredTplCategory(TTemplateCategoryList tplCategoryListIn, CTemplate::ECategory CategoryFilter)
{
    TTemplateCategoryList   filteredTemplateCategoryList;
    for( UInt32 i=0; i<tplCategoryListIn.size(); i++ )
    {
        CTemplate::ECategory category = tplCategoryListIn[i];
        if( category == CTemplate::nCategory_Star )
        {
        }
        else if(CategoryFilter == NSEpic::CTemplate::nCategory_None || CategoryFilter == category)
        {
            filteredTemplateCategoryList.push_back( category );
        }
    }

    return filteredTemplateCategoryList;
}

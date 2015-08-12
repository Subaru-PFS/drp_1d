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

#include <epic/redshift/method/chisquare2solve.h>
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
    m_cut = 1.5;
    m_strongcut = 2.0;

    // Line Matching
    m_minMatchNum = 1;
    m_tol = 0.002;
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
        resultStore.m_dtreepathnum = 1.1;
        {
            //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquare2Solve chiSolve;
            CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
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
        resultStore.m_dtreepathnum = 1.11;
        {
            //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquare2Solve chiSolve;
            CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold );
            if( chisolveResult ) {
                resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
            }
            return true;
        }
    }

    // Ray Match
    Int32 typeFilter = CRay::nType_Emission;
    Int32 restForceFilter = -1;
    if(nRaysDetected == 1){
        restForceFilter =   CRay::nForce_Strong;
    }

    CRayMatching rayMatching;
    CRef<CRayMatchingResult> rayMatchingResult = rayMatching.Compute( rayDetectionResult->RayCatalog, restRayCatalog, redshiftRange, m_minMatchNum, m_tol, typeFilter, -1, restForceFilter);

    if(rayMatchingResult!=NULL){
        // Store matching results
        resultStore.StoreGlobalResult( "raymatching", *rayMatchingResult );

        //check ray matching results
        if(rayMatchingResult->GetSolutionsListOverNumber(0).size()<1){
            Log.LogInfo( "DTreeA - Not match found [1], switching to Chisquare");
            resultStore.m_dtreepathnum = 3.1;
            {  //chisquare
                // Create redshift initial list by spanning redshift acdross the given range, with the given delta
                TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
                DebugAssert( redshifts.size() > 0 );

                CMethodChisquare2Solve chiSolve;
                CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                    tplCatalog, tplCategoryList,
                                                                                    lambdaRange, redshifts, overlapThreshold );
                if( chisolveResult ) {
                    resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
                }
                return true;
            }
        }
    }else{ //TBD: is that likely to happen ?
        Log.LogInfo( "DTreeA - Not match found  [0], switching to Chisquare");
        resultStore.m_dtreepathnum = 3.11;
        {  //chisquare
            // Create redshift initial list by spanning redshift acdross the given range, with the given delta
            TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
            DebugAssert( redshifts.size() > 0 );

            CMethodChisquare2Solve chiSolve;
            CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold );
            if( chisolveResult ) {
                resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
            }
            return true;
        }
    }

    Int32 matchNum = rayMatchingResult->GetMaxMatchingNumber();
    UInt32 nStrongPeaks = rayDetectionResult->RayCatalog.GetFilteredList(CRay::nType_Emission, CRay::nForce_Strong).size();

    //
    if(matchNum >= 1){
        resultStore.m_dtreepathnum = 4.1;
        Log.LogInfo( "DTreeA - compute chisquare on redshift candidates from ray matching" );
        rayMatchingResult->FilterWithRules(spc, lambdaRange, m_winsize);
        TFloat64List redshifts = rayMatchingResult->GetExtendedRedshiftCandidatesOverNumber(0, redshiftStep, 0.01);
        Log.LogInfo( "DTreeA - (n candidates = %d)", redshifts.size());
        { //chisolve with emission templates only
            CMethodChisquare2Solve chiSolve;
            CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, redshifts, overlapThreshold );
            if( chisolveResult ) {
                resultStore.StoreGlobalResult( "redshiftresult", *chisolveResult );
            }
            return true;
        }
    }

    Log.LogInfo( "DTreeA - no other path found than switching to Chisquare...");
    resultStore.m_dtreepathnum = 5.1;
    {  //chisquare
        // Create redshift initial list by spanning redshift acdross the given range, with the given delta
        TFloat64List redshifts = redshiftRange.SpreadOver( redshiftStep );
        DebugAssert( redshifts.size() > 0 );

        CMethodChisquare2Solve chiSolve;
        CConstRef<CChisquare2SolveResult> chisolveResult = chiSolve.Compute( resultStore, spc, spcWithoutCont,
                                                                            tplCatalog, tplCategoryList,
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


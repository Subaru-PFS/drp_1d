#include <epic/redshift/method/dtree7solve.h>

#include <epic/core/log/log.h>
#include <epic/core/debug/assert.h>
#include <epic/redshift/spectrum/template/catalog.h>
#include <epic/redshift/operator/correlation.h>
#include <epic/redshift/operator/chisquare.h>
#include <epic/redshift/extremum/extremum.h>
#include <epic/redshift/processflow/datastore.h>

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
#include <epic/redshift/method/dtree7solve.h>
#include <epic/redshift/method/dtree7solveresult.h>

using namespace NSEpic;
using namespace std;


COperatorDTree7Solve::COperatorDTree7Solve()
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

COperatorDTree7Solve::~COperatorDTree7Solve()
{

}

std::shared_ptr<CDTree7SolveResult> COperatorDTree7Solve::Compute(CDataStore& dataStore, const CSpectrum& spc, const CSpectrum& spcWithoutCont,
                                                        const CTemplateCatalog& tplCatalog, const TStringList& tplCategoryList, const CRayCatalog &restRayCatalog,
                                                        const TFloat64Range& lambdaRange, const TFloat64Range& redshiftRange, Float64 redshiftStep,
                                                        Int32 correlationExtremumCount, Float64 overlapThreshold )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( dataStore, "dtree7solve" );

    dataStore.GetScopedParam( "winsize", m_winsize, 250.0 );
    dataStore.GetScopedParam( "cut", m_cut, 5.0 );
    dataStore.GetScopedParam( "strongcut", m_strongcut, 2.0 );
    dataStore.GetScopedParam( "minMatchNum", m_minMatchNum, 1.0 );
    dataStore.GetScopedParam( "tol", m_tol, 0.002 );


    storeResult = SolveDecisionalTree7(dataStore, spc, spcWithoutCont,
                                       tplCatalog, tplCategoryList, restRayCatalog,
                                       lambdaRange, redshiftRange, redshiftStep,
                                       correlationExtremumCount, overlapThreshold );

    //storeResult = true;
    if( storeResult )
    {
        return std::shared_ptr<CDTree7SolveResult>( new CDTree7SolveResult() );
    }

    return NULL;
}

Bool COperatorDTree7Solve::SolveDecisionalTree7(CDataStore &dataStore, const CSpectrum &spc, const CSpectrum &spcWithoutCont, const CTemplateCatalog &tplCatalog, const TStringList &tplCategoryList, const CRayCatalog &restRayCatalog, const TFloat64Range &lambdaRange, const TFloat64Range &redshiftRange, Float64 redshiftStep, Int32 correlationExtremumCount, Float64 overlapThreshold)
{
    //Log.LogInfo( "Process Decisional Tree 7" );

    //COperatordataStore::CAutoScope resultScope( dataStore, "dtree7" );

    TStringList   filteredTemplateCategoryList = getFilteredTplCategory( tplCategoryList, "emission" );

    CPeakDetection peakDetection(m_winsize, m_cut);
    std::shared_ptr<const CPeakDetectionResult> peakDetectionResult = peakDetection.Compute( spc, lambdaRange);
    dataStore.StoreScopedGlobalResult( "peakdetection", peakDetectionResult );
    Log.LogInfo( "DTree7 - Peak Detection output: %d peaks found", peakDetectionResult->PeakList.size());

    // check Peak Detection results
    if(peakDetectionResult->PeakList.size()<1){
        Log.LogInfo( "DTree7 - No Peak found, switching to Blindsolve");
        m_dtreepathnum = 1.1;
        { //blindsolve
            COperatorBlindSolve blindSolve;
            auto blindsolveResult = blindSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshiftRange, redshiftStep,
                                                                                correlationExtremumCount, overlapThreshold );
            if( blindsolveResult ) {
                dataStore.StoreScopedGlobalResult( "redshiftresult", blindsolveResult );
            }
            return true;
        }
    }

    // Line Detection
    CLineDetection lineDetection(CRay::nType_Emission, m_cut, m_strongcut );
    auto lineDetectionResult = lineDetection.Compute( spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );
    dataStore.StoreScopedGlobalResult( "raycatalog", lineDetectionResult );
    Log.LogInfo( "DTree7 - Line Detection output: %d line(s) found", lineDetectionResult->RayCatalog.GetList().size());

    // check Line Detection results
    Int32 nRaysDetected = lineDetectionResult->RayCatalog.GetList().size();
    if( nRaysDetected < 1){
        Log.LogInfo( "DTree7 - No lines found, switching to ProcessWithoutEL");
        m_dtreepathnum = 1.11;
        { //blindsolve
            COperatorBlindSolve blindSolve;
            auto blindsolveResult = blindSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshiftRange, redshiftStep,
                                                                                correlationExtremumCount, overlapThreshold );
            if( blindsolveResult ) {
                dataStore.StoreScopedGlobalResult( "redshiftresult", blindsolveResult );
            }
            return true;
        }
    }

    // Line Matching
    CRayMatching rayMatching;
    auto rayMatchingResult = rayMatching.Compute(lineDetectionResult->RayCatalog, restRayCatalog, redshiftRange, m_minMatchNum, m_tol );
    if(rayMatchingResult!=NULL){
        // Store matching results
        dataStore.StoreScopedGlobalResult( "raymatching", rayMatchingResult );

        //check line matching results
        if(rayMatchingResult->GetSolutionsListOverNumber(0).size()<1){
            Log.LogInfo( "DTree7 - Not match found [1], switching to ProcessWithoutEL");
            m_dtreepathnum = 1.2;
            { //blindsolve
                COperatorBlindSolve blindSolve;
                auto blindsolveResult = blindSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                                    tplCatalog, tplCategoryList,
                                                                                    lambdaRange, redshiftRange, redshiftStep,
                                                                                    correlationExtremumCount, overlapThreshold );
                if( blindsolveResult ) {
                    dataStore.StoreScopedGlobalResult( "redshiftresult", blindsolveResult );
                }
                return true;
            }
        }
    }else{
        Log.LogInfo( "DTree7 - Not match found  [0], switching to ProcessWithoutEL");
        m_dtreepathnum = 1.21;
        { //blindsolve
            COperatorBlindSolve blindSolve;
            auto blindsolveResult = blindSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                                tplCatalog, tplCategoryList,
                                                                                lambdaRange, redshiftRange, redshiftStep,
                                                                                correlationExtremumCount, overlapThreshold );
            if( blindsolveResult ) {
                dataStore.StoreScopedGlobalResult( "redshiftresult", blindsolveResult );
            }
            return true;
        }
    }

    Int32 matchNum = rayMatchingResult->GetMaxMatchingNumber();
    UInt32 nStrongPeaks = lineDetectionResult->RayCatalog.GetFilteredList(CRay::nType_Emission, CRay::nForce_Strong).size();

    // match num >= 3, or no strong peaks
    if(matchNum >= 3 || (nStrongPeaks < 1 && matchNum >= 2)){
        if(matchNum >= 3){
            Log.LogInfo( "DTree7 - match num >= 3");
            m_dtreepathnum = 2.0;
        }
        if(nStrongPeaks < 1 && matchNum >= 2){
            Log.LogInfo( "DTree7 - n Strong Peaks < 1, MatchNum>=2");
            m_dtreepathnum = 2.1;
        }
        Log.LogInfo( "DTree7 - compute merits on redshift candidates from line matching" );
        TFloat64List roundedRedshift = rayMatchingResult->GetRoundedRedshiftCandidatesOverNumber(matchNum-1, redshiftStep);
        Log.LogInfo( "DTree7 - (n candidates = %d)", roundedRedshift.size());
        { //chisolve with emission templqtes only
            CMethodChisquareSolve chiSolve;
            auto chisolveResult = chiSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, roundedRedshift, overlapThreshold );
            if( chisolveResult ) {
                dataStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
            }
            return true;
        }
    }

    //    // 3 lines 1 match or 4 lines 2 matches
    //    if(nRaysDetected - matchNum >= 2){
    //        Log.LogInfo( "3 lines 1 match or 4 lines 2 matches...");
    //    }

    // use Strong peaks
    if(nStrongPeaks > 0){
        Log.LogInfo( "DTree7 - Line Matching with %d strong peaks", nStrongPeaks);
        CRayMatching rayMatchingStrong;
        auto rayMatchingStrongResult = rayMatchingStrong.Compute(lineDetectionResult->RayCatalog, restRayCatalog, redshiftRange, m_minMatchNum, m_tol, CRay::nType_Emission, CRay::nForce_Strong );
        Int32 matchNumStrong = rayMatchingStrongResult->GetMaxMatchingNumber();

        if(matchNumStrong>1){
            Log.LogInfo( "DTree7 - match num strong >= 2, compute merits on redshift candidates from strong line matching");
            TFloat64List roundedRedshift = rayMatchingStrongResult->GetRoundedRedshiftCandidatesOverNumber(matchNumStrong-1, redshiftStep);
            m_dtreepathnum = 2.2;
            { //chisolve with emission templqtes only
                CMethodChisquareSolve chiSolve;
                auto chisolveResult = chiSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                                    tplCatalog, filteredTemplateCategoryList,
                                                                                    lambdaRange, roundedRedshift, overlapThreshold );
                if( chisolveResult ) {
                    dataStore.StoreScopedGlobalResult( "redshiftresult", chisolveResult );
                }
                return true;
            }
        }else{
            Log.LogInfo( "DTree7 - Not match found with strong lines, switching to ProcessWithoutEL (EZ: only_correlation_... equivalent)");
            m_dtreepathnum = 3.0;
            { // corrsolve with emission templates only
                COperatorCorrelationSolve Solve;
                auto solveResult = Solve.Compute( dataStore, spc, spcWithoutCont,
                                                                                tplCatalog, filteredTemplateCategoryList,
                                                                                lambdaRange, redshiftRange, redshiftStep,
                                                                                overlapThreshold );
                if( solveResult ) {
                    dataStore.StoreScopedGlobalResult( "redshiftresult", solveResult );
                }
                return true;
            }
        }
    }

    Log.LogInfo( "DTree7 - no other path found than switching to ProcessWithoutEL...");
    m_dtreepathnum = 1.3;
    { //blindsolve
        COperatorBlindSolve blindSolve;
        auto blindsolveResult = blindSolve.Compute( dataStore, spc, spcWithoutCont,
                                                                            tplCatalog, tplCategoryList,
                                                                            lambdaRange, redshiftRange, redshiftStep,
                                                                            correlationExtremumCount, overlapThreshold );
        if( blindsolveResult ) {
            dataStore.StoreScopedGlobalResult( "redshiftresult", blindsolveResult );
        }
        return true;
    }
}


TStringList COperatorDTree7Solve::getFilteredTplCategory( const TStringList& tplCategoryListIn, const std::string& CategoryFilter)
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

#include <RedshiftLibrary/method/linematchingsolve.h>

#include <RedshiftLibrary/debug/assert.h>
#include <RedshiftLibrary/spectrum/template/catalog.h>
#include <RedshiftLibrary/operator/correlation.h>
#include <RedshiftLibrary/operator/chisquare.h>
#include <RedshiftLibrary/extremum/extremum.h>
#include <RedshiftLibrary/processflow/datastore.h>

#include <RedshiftLibrary/operator/peakdetection.h>
#include <RedshiftLibrary/operator/raydetection.h>
#include <RedshiftLibrary/operator/raymatching.h>

#include <RedshiftLibrary/operator/peakdetectionresult.h>
#include <RedshiftLibrary/operator/raydetectionresult.h>
#include <RedshiftLibrary/operator/raymatchingresult.h>
#include <RedshiftLibrary/method/linematchingsolveresult.h>

using namespace NSEpic;
using namespace std;


COperatorLineMatchingSolve::COperatorLineMatchingSolve()
{
//     Peak Detection
//    m_winsize = 250.0;
//    m_cut = 5.0;
//    m_strongcut = 2.0;
//    m_minsize = 3;
//    m_maxsize = 70;
//    m_enlargeRate = 2.0;


//     Line Matching
//    m_minMatchNum = 1;
//    m_tol = 0.002;

//    //eudlid overrides
//    if(0)
//    {
//        m_winsize = 250.0;
//        m_cut = 3.0;
//        m_maxsize = 120;
//    }

//    //pfs TF overrides
//    if(1)
//    {
//        m_winsize = 250.0;
//        m_cut = 150;
//        m_maxsize = 120;
//        m_enlargeRate = 1.0;
//        m_tol = 0.0025;
//    }
}

COperatorLineMatchingSolve::~COperatorLineMatchingSolve()
{

}


const std::string COperatorLineMatchingSolve::GetDescription()
{
    std::string desc;

    desc = "Method linematching:\n";

    desc.append("\tparam: peakdetection.winsize = <float value in Angstrom>\n");
    desc.append("\tparam: peakdetection.cut = <float value>\n");
    desc.append("\tparam: peakdetection.enlargerate = <float value>\n");

    desc.append("\tparam: linedetection.cut = <float value>\n");
    desc.append("\tparam: linedetection.strongcutfactor = <float value>\n");
    desc.append("\tparam: linedetection.minlinewidth = <float value in Angstrom>\n");
    desc.append("\tparam: linedetection.maxlinewidth = <float value in Angstrom>\n");

    desc.append("\tparam: linematching.tolerance = <float value>\n");
    desc.append("\tparam: linematching.minmatchnum = <int value>\n");


    return desc;

}

std::shared_ptr<const CLineMatchingSolveResult> COperatorLineMatchingSolve::Compute(  CDataStore& dataStore, const CSpectrum& spc,
                                                                      const TFloat64Range& lambdaRange, const TFloat64Range& redshiftsRange, Float64 redshiftStep, const CRayCatalog& restRayCatalog )
{
    Bool storeResult = false;

    CDataStore::CAutoScope resultScope( dataStore, "linematchingsolve" );

    Float64 opt_winsize;
    dataStore.GetScopedParam( "peakdetection.winsize", opt_winsize, 250.0 );
    Float64 opt_peakdetectioncut;
    dataStore.GetScopedParam( "peakdetection.cut", opt_peakdetectioncut, 5.0 );
    Float64 opt_peakenlargerate;
    dataStore.GetScopedParam( "peakdetection.enlargerate", opt_peakenlargerate, 2.0 );

    CPeakDetection peakDetection(opt_winsize, opt_peakdetectioncut, 1, opt_peakenlargerate);
    auto peakDetectionResult = peakDetection.Compute( spc, lambdaRange);
    if( peakDetectionResult )
    {
        dataStore.StoreScopedGlobalResult( "peakdetection", peakDetectionResult );
    }else
    {
      return NULL;
    }

    Float64 opt_linedetectioncut;
    dataStore.GetScopedParam( "linedetection.cut", opt_linedetectioncut, 5.0 );
    Float64 opt_strongcutfactor;
    dataStore.GetScopedParam( "linedetection.strongcutfactor", opt_strongcutfactor, 2.0 );
    Float64 opt_minlinewidth;
    dataStore.GetScopedParam( "linedetection.minlinewidth", opt_minlinewidth, 3.0 );
    Float64 opt_maxlinewidth;
    dataStore.GetScopedParam( "linedetection.maxlinewidth", opt_maxlinewidth, 70.0 );

    CLineDetection lineDetection(CRay::nType_Emission, opt_linedetectioncut, opt_strongcutfactor, opt_winsize, opt_minlinewidth, opt_maxlinewidth);
    auto lineDetectionResult = lineDetection.Compute( spc, lambdaRange, peakDetectionResult->PeakList, peakDetectionResult->EnlargedPeakList );

    if( lineDetectionResult ) {
        dataStore.StoreScopedGlobalResult( "raycatalog", lineDetectionResult );

        if(lineDetectionResult->RayCatalog.GetList().size()<1){
            //return NULL;
        }
    }else
    {
      return NULL;
    }

    // --- Match
    Float64 opt_tolerance;
    dataStore.GetScopedParam( "linematching.tolerance", opt_tolerance, 0.002 );
    Float64 opt_minmatchnum;
    dataStore.GetScopedParam( "linematching.minmatchnum", opt_minmatchnum, 1.0 );

    CRayMatching rayMatching;
    auto rayMatchingResult = rayMatching.Compute(lineDetectionResult->RayCatalog, restRayCatalog, redshiftsRange, opt_minmatchnum, opt_tolerance );

    // Store matching results
    if( rayMatchingResult )
        dataStore.StoreScopedGlobalResult( "raymatching", rayMatchingResult );


    storeResult = true; //always save a matching result
    if( storeResult )
    {
        return std::shared_ptr<const CLineMatchingSolveResult>( new CLineMatchingSolveResult() );
    }

    return NULL;
}

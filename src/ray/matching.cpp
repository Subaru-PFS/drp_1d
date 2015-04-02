#include <epic/redshift/ray/matching.h>

#include <algorithm>    // std::sort
#include <math.h>
#include <stdio.h>

using namespace NSEpic;
IMPLEMENT_MANAGED_OBJECT(CRayMatching)

// This class's initial implementation duplicates the VIPGI method to match lines found with catalog lines, as in EZELmatch.py, Line 264

CRayMatching::CRayMatching()
{
}

CRayMatching::~CRayMatching()
{
}

Bool CRayMatching::Compute( const CRayCatalog::TRayVector& detectedRayList, const CRayCatalog::TRayVector& restRayList, const TFloat64Range& redshiftRange, Int32 nThreshold, Float64 tol )
{
    TRedshiftSolutionSetList solutions;
    Int32 nDetectedRay = detectedRayList.size();
    if( nDetectedRay == 1 )
        // if there is only one detected line, then there are N=#restlines possible redshifts
        for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
        {
            TRedshiftSolutionSet solution;
            Float64 redShift=(detectedRayList[0].GetPosition()-restRayList[0].GetPosition())/restRayList[0].GetPosition();
            if( redShift > 0 ){
                solution.push_back(SRedshiftSolution( detectedRayList[0].GetPosition(), restRayList[0].GetPosition(), redShift));
                solutions.push_back(solution);
            }
        }
    else{
        for( UInt32 iDetectedRay=0; iDetectedRay<detectedRayList.size(); iDetectedRay++ )
        {
            for( UInt32 iRestRay=0; iRestRay<restRayList.size(); iRestRay++ )
            {
                // for each detected line / rest line couple, enumerate how many other couples fit within the tolerance
                TRedshiftSolutionSet solution;
                Float64 redShift=(detectedRayList[iDetectedRay].GetPosition()-restRayList[iRestRay].GetPosition())/restRayList[iRestRay].GetPosition();
                if( redShift < 0 ){
                    continue;
                }
                solution.push_back( SRedshiftSolution( detectedRayList[iDetectedRay].GetPosition(), restRayList[iRestRay].GetPosition(), redShift) );

                for( UInt32 iDetectedRay2=0; iDetectedRay2<detectedRayList.size(); iDetectedRay2++ )
                {
                    if(iDetectedRay!=iDetectedRay2){
                        for( UInt32 iRestRay2=0; iRestRay2<restRayList.size(); iRestRay2++ )
                        {
                            Float64 redShift2=(detectedRayList[iDetectedRay2].GetPosition()-restRayList[iRestRay2].GetPosition())/restRayList[iRestRay2].GetPosition();
                            if( redShift2 < 0 ){
                                continue;
                            }
                            Float64 redshiftTolerance = tol*(1+(redShift+redShift2)*0.5);
                            if( fabs( (redShift-redShift2) ) <= redshiftTolerance){
                                Bool found = false;
                                //avoid repeated solution sets
                                for( UInt32 iSet=0; iSet<solution.size(); iSet++ )
                                {
                                    if(solution[iSet].DetectedRay == detectedRayList[iDetectedRay2].GetPosition()){
                                        found = true;
                                        break;
                                    }
                                }
                                if(!found){
                                    solution.push_back( SRedshiftSolution( detectedRayList[iDetectedRay2].GetPosition(), restRayList[iRestRay2].GetPosition(), redShift2) );
                                }
                            }
                        }
                    }
                }
                if(solution.size()>0){
                    solutions.push_back(solution);
                }
            }
        }
    }

    //delete duplicated solutions
    TRedshiftSolutionSetList newSolutions;
    for( UInt32 iSol=0; iSol<solutions.size(); iSol++ )
    {
        TRedshiftSolutionSet currentSet = solutions[iSol];
        sort(currentSet.begin(), currentSet.end());
        Bool found = false;
        for( UInt32 iNewSol=0; iNewSol<newSolutions.size(); iNewSol++ )
        {
            if( AreSolutionSetsEqual(newSolutions[iNewSol],currentSet) ){
                found = true;
            }
        }
        if(!found)
        {
            if(currentSet.size() > nThreshold){
                Float64 redshiftMean = 0.f;
                for( UInt32 i=0; i<currentSet.size(); i++ )
                {
                    redshiftMean += currentSet[i].Redshift;
                }
                redshiftMean /= currentSet.size();
                if(redshiftMean>redshiftRange.GetBegin() && redshiftMean<redshiftRange.GetEnd()){
                    newSolutions.push_back(currentSet);
                }
            }
        }
    }

    return true;
}

Bool CRayMatching::AreSolutionSetsEqual(TRedshiftSolutionSet s1, TRedshiftSolutionSet s2)
{
    if(s1.size() != s2.size()){
        return false;
    }
    Bool diffFound = false;
    for( UInt32 iSet=0; iSet<s1.size(); iSet++ )
    {
        if(s1[iSet].DetectedRay!=s2[iSet].DetectedRay || s1[iSet].RestRay!=s2[iSet].RestRay || s1[iSet].Redshift!=s2[iSet].Redshift){
           diffFound = true;
           return false;
        }
    }

    if(!diffFound){
        return true;
    }
}

const Float64 CRayMatching::GetResults() const
{
    return m_Results;
}



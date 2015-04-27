#include <epic/redshift/operator/raymatchingresult.h>

using namespace NSEpic;

IMPLEMENT_MANAGED_OBJECT( CRayMatchingResult )

CRayMatchingResult::CRayMatchingResult()
{

}

CRayMatchingResult::~CRayMatchingResult()
{

}

Void CRayMatchingResult::Save( std::ostream& stream ) const
{

}


Bool CRayMatchingResult::GetBestRedshift(Float64& Redshift, Int32& MatchingNumber) const
{
    MatchingNumber = GetMaxMatchingNumber();
    TSolutionSetList selectedResults = GetSolutionsListOverNumber(MatchingNumber-1);

    if(selectedResults.size()>0){
        Redshift = GetMeanRedshiftSolution(selectedResults[0]);
    }else{
        return false;
    }

    return true;
}


CRayMatchingResult::TSolutionSetList CRayMatchingResult::GetSolutionsListOverNumber(Int32 number) const
{
    //select results by matching number
    TSolutionSetList selectedResults;
    for( UInt32 iSol=0; iSol<SolutionSetList.size(); iSol++ )
    {
        TSolutionSet currentSet = SolutionSetList[iSol];
        if(currentSet.size() > number){
            selectedResults.push_back(currentSet);
        }

    }
    return selectedResults;
}

Float64 CRayMatchingResult::GetMeanRedshiftSolutionByIndex(Int32 index) const
{
    if(index > SolutionSetList.size()-1)
    {
        return -1.0;
    }

    Float64 redshiftMean = 0.0;
    const TSolutionSet& currentSet = SolutionSetList[index];
    for( UInt32 i=0; i<currentSet.size(); i++ )
    {
        redshiftMean += currentSet[i].Redshift;
    }
    redshiftMean /= currentSet.size();

    return redshiftMean;
}


Float64 CRayMatchingResult::GetMeanRedshiftSolution( const TSolutionSet& s ) const
{
    TSolutionSet currentSet = s;
    Float64 redshiftMean=0.0;
    for( UInt32 i=0; i<currentSet.size(); i++ )
    {
        redshiftMean += currentSet[i].Redshift;
    }
    redshiftMean /= (float)currentSet.size();

    return redshiftMean;
}

Int32 CRayMatchingResult::GetMaxMatchingNumber() const
{
    if(SolutionSetList.size() < 1)
    {
        return -1.0;
    }

    Int32 maxNumber = 0;
    for( UInt32 i=0; i<SolutionSetList.size() ; i++ )
    {
        TSolutionSet currentSet = SolutionSetList[i];
        if(maxNumber<currentSet.size()){
            maxNumber = currentSet.size();
        }

    }

    return maxNumber;
}

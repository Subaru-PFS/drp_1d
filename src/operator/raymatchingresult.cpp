#include <epic/redshift/operator/raymatchingresult.h>
#include <stdio.h>

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
    TSolutionSetList selectedResults = GetSolutionsListOverNumber(0);

    if(selectedResults.size()>0){
        std::string strList;
        char tmpChar[256];
        strList.append("#MATCH_NUM\tDETECTED_LINES\tREST_LINES\tZ\n");
        for( UInt32 iSol=0; iSol<selectedResults.size(); iSol++ )
        {
            std::string strMNUM= "";
            TSolutionSet currentSet = selectedResults[iSol];
            sprintf(tmpChar, "%d", currentSet.size());
            strMNUM.append(tmpChar);
            std::string strZ= "";
            sprintf(tmpChar, "%.10f", GetMeanRedshiftSolution(currentSet) );
            strZ.append(tmpChar);

            std::string strDetected = "(";
            std::string strRest = "(";
            for( UInt32 i=0; i<currentSet.size(); i++ )
            {
                sprintf(tmpChar, "%.1f, ", currentSet[i].DetectedRay);
                strDetected.append(tmpChar);
                sprintf(tmpChar, "%.1f, ", currentSet[i].RestRay);
                strRest.append(tmpChar);
            }
            strDetected = strDetected.substr(0, strDetected.size()-2);
            strDetected.append(")");
            strRest = strRest.substr(0, strRest.size()-2);
            strRest.append(")");

            strList.append(strMNUM);
            strList.append("\t");
            strList.append(strDetected);
            strList.append("\t");
            strList.append(strRest);
            strList.append("\t");
            strList.append(strZ);
            strList.append("\n");
        }
        stream << strList.c_str();
    }

    return;
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

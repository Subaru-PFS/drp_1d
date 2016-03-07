#include <cstdarg>
#include <iostream>

#include <epic/core/log/log.h>
#include <epic/redshift/ray/ruleStrongHigherThanWeak.h>

using namespace NSEpic;
using namespace std;

void CRuleStrongHigherThanWeak::SetUp( Bool EnabledArgument, ... )
{
  Name = "strongweak";
  Enabled = EnabledArgument;
  va_list Arguments;
  va_start ( Arguments, EnabledArgument );
  m_LineType = va_arg( Arguments, Int32 );
}

/** \brief Verify that "stronger lines have higher amplitudes than weaker lines" rule is applicable, and then apply it.
 * If the maximum amplitude for strong lines of the specified type is -1, do nothing.
 * For each weak line of the specified type:
 *  Find the first element that contains that weak line
 *  Get the index of the entry in that element that corresponds to the weak line
 *  If the indexed entry IsOutsideLambdaRange, go for the next weak line
 *  Get the parameters for the entry
 *  Limit the amplitude of the entry to the maximum amplitude for strong lines
 **/
void CRuleStrongHigherThanWeak::Correct( CLineModelElementList& LineModelElementList )
{
  Float64 coeff = 1.0;
  Float64 erStrong=-1.0;
  Float64 maxiStrong = FindHighestStrongLineAmp( m_LineType, erStrong, LineModelElementList );
  if( maxiStrong == -1 )
    {
      //Log.LogDebug( "Rule %s: no strong line detected.", Name.c_str() );
      return;
    }
  for( UInt32 iRestRayWeak=0; iRestRayWeak<LineModelElementList.m_RestRayList.size(); iRestRayWeak++ ) //loop on the weak lines
    {
      if( LineModelElementList.m_RestRayList[iRestRayWeak].GetForce() != CRay::nForce_Weak )
	{
	  continue;
        }
      if( LineModelElementList.m_RestRayList[iRestRayWeak].GetType() != m_LineType )
	{
	  continue;
	}
      Int32 eIdxWeak = LineModelElementList.FindElementIndex( iRestRayWeak );
      Int32 subeIdxWeak = LineModelElementList.m_Elements[eIdxWeak]->FindElementIndex( iRestRayWeak );
      if( LineModelElementList.m_Elements[eIdxWeak]->IsOutsideLambdaRange( subeIdxWeak ) == true )
	{
	  continue;
	}
      //Log.LogDebug( "Rule %s: element %d has force weak, type %d and is not outside lambda range.", Name.c_str(), iRestRayWeak, m_LineType );
      Float64 nSigma = 1.0;
      Float64 ampA = maxiStrong;
      Float64 erA = erStrong;
      Float64 ampB = LineModelElementList.m_Elements[eIdxWeak]->GetFittedAmplitude( subeIdxWeak );
      Float64 erB = LineModelElementList.m_Elements[eIdxWeak]->GetFittedAmplitudeErrorSigma( subeIdxWeak );
      Float64 maxB = (coeff*ampA) + coeff*(erA*nSigma);
      LineModelElementList.m_Elements[eIdxWeak]->LimitFittedAmplitude( subeIdxWeak, maxB );
    }
}

Bool CRuleStrongHigherThanWeak::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

/**
 * \brief Returns the maximum amplitude between strong lines within the support of m_Elements. The referenced er argument will hold the error sigma for the same element.
 **/
Float64 CRuleStrongHigherThanWeak::FindHighestStrongLineAmp( Int32 linetype , Float64 &er, CLineModelElementList& LineModelElementList )
{
  Float64 maxi = -1.0;
  for( UInt32 iRestRayStrong=0; iRestRayStrong<LineModelElementList.m_RestRayList.size(); iRestRayStrong++ ) //loop on the strong lines
    {
      if( LineModelElementList.m_RestRayList[iRestRayStrong].GetForce() != CRay::nForce_Strong )
	{
	  continue;
	}
      if( LineModelElementList.m_RestRayList[iRestRayStrong].GetType() != m_LineType )
	{
	  continue;
	}
      Int32 eIdxStrong = LineModelElementList.FindElementIndex( iRestRayStrong );
      Int32 subeIdxStrong = LineModelElementList.m_Elements[eIdxStrong]->FindElementIndex( iRestRayStrong );
      if( LineModelElementList.m_Elements[eIdxStrong]->IsOutsideLambdaRange( subeIdxStrong ) == true )
	{
	  continue;
	}
      Float64 ampStrong = LineModelElementList.m_Elements[eIdxStrong]->GetFittedAmplitude( subeIdxStrong );
      if( maxi<ampStrong )
	{
	  maxi = ampStrong;
	  er = LineModelElementList.m_Elements[eIdxStrong]->GetFittedAmplitudeErrorSigma( subeIdxStrong );
	}
    }
  //Log.LogDebug( "Highest strong line amplitude = %f", maxi );
  return maxi;
}

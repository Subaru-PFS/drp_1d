#include <cstdarg>
#include <iostream>

#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/ruleStrongHigherThanWeak.h>

using namespace NSEpic;
using namespace std;

CRuleStrongHigherThanWeak::CRuleStrongHigherThanWeak():
  m_LineType(0)
{
}

CRuleStrongHigherThanWeak::~CRuleStrongHigherThanWeak()
{
}

void CRuleStrongHigherThanWeak::SetUp( Bool EnabledArgument, ... )
{
  Name = "strongweak";
  Enabled = EnabledArgument;
  va_list Arguments;
  va_start ( Arguments, EnabledArgument );
  m_LineType = va_arg( Arguments, Int32 );
  va_end(Arguments);
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
  std::string strongName = "";
  Float64 maxiStrong = FindHighestStrongLineAmp( m_LineType, erStrong, strongName, LineModelElementList );
  if( maxiStrong == -1 )
    {
      //Log.LogDebug( "Rule %s: no strong line detected.", Name.c_str() );

      //case 0
      return;

      //case 1:
      //maxiStrong = 0.0;
      //erStrong = 0.0;
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

      //Method 0 : no noise taken into acccount
      //Float64 maxB = (coeff*ampA);

      //Method 1 : using Strong line noise to be more tolerant
      Float64 maxB = (coeff*ampA) + coeff*(erA*nSigma);

      //Method 2 : using Strong line noise and Weak line noise to correct with a ratio
//      Float64 maxB = ampB; //default value
//      if(erB>0.0 && erB<erA && erA>0.0)
//      {
//          maxB = (coeff*ampA)*(erA/erB);
//      }else{
//          maxB = (coeff*ampA);
//      }
      //

      if(maxB==std::min(maxB, ampB) && maxB!=ampB)
      {
          LineModelElementList.m_Elements[eIdxWeak]->LimitFittedAmplitude( subeIdxWeak, maxB );
          //log the correction
          {
              std::string nameWeak = LineModelElementList.m_RestRayList[iRestRayWeak].GetName();
              if(Logs.size()==0)
              {
                  std::string strTmp0 = boost::str( (boost::format("correct - %-10s") % "STRONG_WEAK" ));
                  Logs.append(strTmp0.c_str());
              }
              std::string strTmp = boost::str( (boost::format("\n\tlineWeak=%-10s, lineStrong=%-10s, previousAmp=%.4e, correctedAmp=%.4e") % nameWeak % strongName % ampB % maxB ));
              Logs.append(strTmp.c_str());
          }
      }
    }
}

Bool CRuleStrongHigherThanWeak::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

/**
 * \brief Returns the maximum amplitude between superstrong lines within the support of m_Elements. The referenced er argument will hold the error sigma for the same element.
 **/
Float64 CRuleStrongHigherThanWeak::FindHighestStrongLineAmp( Int32 linetype , Float64 &er, std::string &name, CLineModelElementList& LineModelElementList )
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
      Float64 erStrong = LineModelElementList.m_Elements[eIdxStrong]->GetFittedAmplitudeErrorSigma( subeIdxStrong );
      // if(erStrong>0.0 && ampStrong>0.0)
      // {
      //     lineSnr = ampStrong/erStrong;
      // }
      if( maxi<ampStrong /*&& lineSnr>validSNRCut*/ )
	{
	  maxi = ampStrong;
      er = erStrong;
      name = LineModelElementList.m_RestRayList[iRestRayStrong].GetName();
    }
    }
  //Log.LogDebug( "Highest strong line amplitude = %f", maxi );
  return maxi;
}

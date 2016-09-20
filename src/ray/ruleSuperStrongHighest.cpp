#include <cstdarg>
#include <iostream>

#include <epic/core/log/log.h>
#include <epic/redshift/ray/ruleSuperStrongHighest.h>

using namespace NSEpic;
using namespace std;

void CRuleSuperStrong::SetUp( Bool EnabledArgument, ... )
{
  Name = "superstrong";
  Enabled = EnabledArgument;
  va_list Arguments;
  va_start ( Arguments, EnabledArgument );
  m_LineType = va_arg( Arguments, Int32 );
  std::string _superStrongTag1 = std::string ( va_arg( Arguments, const char* ) );
  m_SuperStrongTags.push_back(_superStrongTag1);
  std::string _superStrongTag2 = std::string ( va_arg( Arguments, const char* ) );
  m_SuperStrongTags.push_back(_superStrongTag2);
  std::string _superStrongTag3 = std::string ( va_arg( Arguments, const char* ) );
  m_SuperStrongTags.push_back(_superStrongTag3);
}

/** \brief Verify that "super stronger lines have higher amplitudes than other lines" rule is applicable, and then apply it.
 * If the maximum amplitude for super strong lines of the specified type is -1, do nothing.
 * For each other line of the specified type:
 *  Find the first element that contains that line
 *  Get the index of the entry in that element that corresponds to the line
 *  If the indexed entry IsOutsideLambdaRange, go for the next line
 *  Get the parameters for the entry
 *  Limit the amplitude of the entry to the maximum amplitude for super strong lines
 **/
void CRuleSuperStrong::Correct( CLineModelElementList& LineModelElementList )
{
  Float64 coeff = 1.0;
  Float64 erStrong=-1.0;
  std::string strongName = "";
  Float64 maxiStrong = FindHighestSuperStrongLineAmp( m_SuperStrongTags, erStrong, strongName, LineModelElementList );
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

      bool foundSuperStringTag=false;
      for(Int32 k=0; k<m_SuperStrongTags.size(); k++)
      {
          if(LineModelElementList.m_RestRayList[iRestRayWeak].GetName() == m_SuperStrongTags[k])
          {
              foundSuperStringTag=true;
              break;
          }
      }
      if( foundSuperStringTag )
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
                  std::string strTmp0 = boost::str( (boost::format("correct - %-10s") % "SUPER_STRONG" ));
                  Logs.append(strTmp0.c_str());
              }
              std::string strTmp = boost::str( (boost::format("\n\tline=%-10s, lineSuperStrong=%-10s, previousAmp=%.4e, correctedAmp=%.4e") % nameWeak % strongName % ampB % maxB ));
              Logs.append(strTmp.c_str());
          }
      }
    }
}

Bool CRuleSuperStrong::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

/**
 * \brief Returns the maximum amplitude between strong lines within the support of m_Elements. The referenced er argument will hold the error sigma for the same element.
 **/
Float64 CRuleSuperStrong::FindHighestSuperStrongLineAmp( TStringList superstrongTags , Float64 &er, std::string &name, CLineModelElementList& LineModelElementList )
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
      bool foundSuperStringTag=false;
      for(Int32 k=0; k<superstrongTags.size(); k++)
      {
        if(LineModelElementList.m_RestRayList[iRestRayStrong].GetName() == superstrongTags[k])
        {
            foundSuperStringTag=true;
            break;
        }
      }
      if( !foundSuperStringTag )
      {
        continue;
      }

      Int32 eIdxStrong = LineModelElementList.FindElementIndex( iRestRayStrong );
      Int32 subeIdxStrong = LineModelElementList.m_Elements[eIdxStrong]->FindElementIndex( iRestRayStrong );
      if( LineModelElementList.m_Elements[eIdxStrong]->IsOutsideLambdaRange( subeIdxStrong ) == true )
	{
	  continue;
	}

      Float64 validSNRCut = 0.05;
      Float64 ampStrong = LineModelElementList.m_Elements[eIdxStrong]->GetFittedAmplitude( subeIdxStrong );
      Float64 erStrong = LineModelElementList.m_Elements[eIdxStrong]->GetFittedAmplitudeErrorSigma( subeIdxStrong );
      Float64 lineSnr = -1.0;
      if(erStrong>0.0 && ampStrong>0.0)
      {
          lineSnr = ampStrong/erStrong;
      }
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

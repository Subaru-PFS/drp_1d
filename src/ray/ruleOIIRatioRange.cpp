#include <cstdarg>
#include <iostream>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <epic/core/log/log.h>
#include <epic/redshift/ray/ruleOIIRatioRange.h>

using namespace NSEpic;
using namespace std;

void CRuleOIIRatioRange::SetUp( Bool EnabledArgument, ... )
{
  Name = "oiiratio";
  Enabled = EnabledArgument;
  va_list Arguments;
  va_start ( Arguments, EnabledArgument );
  m_LineType = va_arg( Arguments, Int32 );
  m_LineA = std::string ( va_arg( Arguments, const char* ) );
  m_LineB = std::string ( va_arg( Arguments, const char* ) );
  m_Coefficient = va_arg( Arguments, Float64 );
}

Bool CRuleOIIRatioRange::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

/**
 * For two distinct lines, if neither IsOutsideLambdaRange, and their amplitudes are beyond a range (considering coeff), SetFittedAmplitude of each with corrected values.
 **/
Void CRuleOIIRatioRange::Correct( CLineModelElementList& LineModelElementList )
{
  Int32 iA = LineModelElementList.FindElementIndex( m_LineA, m_LineType);
  if( iA==-1 )
    {
      Log.LogDebug( "Rule %s: line %s not found.", Name.c_str(), m_LineA.c_str() );
      return;
    }
  if( LineModelElementList.m_Elements[iA]->GetSize()>1 )
    {
      Log.LogDebug( "Rule %s: line %s has size < 1.", Name.c_str(), m_LineA.c_str() );
      iA=-1;
    }
  Int32 iB = LineModelElementList.FindElementIndex( m_LineB, m_LineType );
  if( iB==-1 )
    {
      Log.LogDebug( "Rule %s: line %s not found.", Name.c_str(), m_LineB.c_str() );
      return;
    }
  if( LineModelElementList.m_Elements[iB]->GetSize()>1 )
    {
      Log.LogDebug( "Rule %s: line %s has size < 1.", Name.c_str(), m_LineB.c_str() );
      iB=-1;
    }
  if (iA==-1 || iB==-1 || iA==iB )
    {
      return;
    }
  if( LineModelElementList.m_Elements[iA]->IsOutsideLambdaRange() == false && LineModelElementList.m_Elements[iB]->IsOutsideLambdaRange() == false)
    {
      Float64 ampA = LineModelElementList.m_Elements[iA]->GetFittedAmplitude( 0 );
      Float64 erA = LineModelElementList.m_Elements[iA]->GetFittedAmplitudeErrorSigma( 0 );
      Float64 ampB = LineModelElementList.m_Elements[iB]->GetFittedAmplitude( 0 );
      Float64 erB = LineModelElementList.m_Elements[iB]->GetFittedAmplitudeErrorSigma( 0 );
      Int32 i1 = iA;
      Int32 i2 = iB;
      Float64 amp1 = ampA;
      Float64 er1 = erA;
      Float64 amp2 = ampB;
      Float64 er2 = erB;
      if( std::abs( ampA ) > std::abs( ampB*m_Coefficient ) )
	{
	  i1 = iA;
	  i2 = iB;
	  amp1 = ampA;
	  er1 = erA;
	  amp2 = ampB;
	  er2 = erB;
        }
      else
	{
	  if( std::abs(ampB) > std::abs(ampA*m_Coefficient) )
	    {
	      i1 = iB;
	      i2 = iA;
	      amp1 = ampB;
	      er1 = erB;
	      amp2 = ampA;
	      er2 = erA;
	    }
	  else
	    {
            return;
	    }
	}
      Float64 R = m_Coefficient;
      Float64 w1 = 0.0;
      if( er1!=0.0 )
	{
	  w1 = 1.0/(er1*er1);
        }
      Float64 w2 = 0.0;
      if( er2!=0.0 )
	{
	  w2 = 1.0/(er2*er2*R*R);
        }
      Float64 corrected1 = (amp1*w1 + amp2*w2*R)/(w1+w2);
      Float64 corrected2 = corrected1/R;

      {
          //log the correction
          Float64 correctedA = corrected1;
          Float64 correctedB = corrected2;
          if(i1==iB)
          {
              correctedA = corrected2;
              correctedB = corrected1;
          }
          std::string strTmp0 = boost::str( (boost::format("correct - %-10s") % "RATIO_RANGE" ));
          Logs.append(strTmp0.c_str());
          std::string strTmp1 = boost::str( (boost::format("\n\tline=%-10s, previousAmp=%.4e, correctedAmp=%.4e") % m_LineA % ampA % correctedA) );
          Logs.append(strTmp1.c_str());
          std::string strTmp2 = boost::str( (boost::format("\n\tline=%-10s, previousAmp=%.4e, correctedAmp=%.4e") % m_LineB % ampB % correctedB) );
          Logs.append(strTmp2.c_str());
      }

      LineModelElementList.m_Elements[i1]->SetFittedAmplitude( corrected1, er1 );
      LineModelElementList.m_Elements[i2]->SetFittedAmplitude( corrected2, er2 );
    }
}

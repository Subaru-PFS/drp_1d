// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
// 
// https://www.lam.fr/
// 
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
// 
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info". 
// 
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 
// 
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 
// 
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#include <cstdarg>
#include <iostream>

#include "RedshiftLibrary/log/log.h"
#include "RedshiftLibrary/ray/rule2singlelinesamplitude.h"
#include "boost/format.hpp"

using namespace NSEpic;
using namespace std;

CRule2SingleLinesAmplitude::CRule2SingleLinesAmplitude() :
  m_LineType(0),
  m_LineA(""),
  m_LineB(""),
  m_Coefficient(0)
{
}


void CRule2SingleLinesAmplitude::SetUp( bool EnabledArgument, ... )
{
  Name = "balmersingle";
  Enabled = EnabledArgument;
  va_list Arguments;
  va_start ( Arguments, EnabledArgument );
  m_LineType = va_arg( Arguments, Int32 );
  m_LineA = std::string ( va_arg( Arguments, const char* ) );
  m_LineB = std::string ( va_arg( Arguments, const char* ) );
  m_Coefficient = va_arg( Arguments, Float64 );

  if(0){
      std::string strTmp = boost::str( (boost::format("setup - linetype=%d, lineA=%s, lineB=%s") % m_LineType % m_LineA % m_LineB) );
      Logs.append(strTmp.c_str());
  }
  va_end(Arguments);
}

/**
 * \brief Correct both lines depending on their sigmas.
 **/
void CRule2SingleLinesAmplitude::Correct( CLineModelElementList& LineModelElementList )
{
    Int32 iA = LineModelElementList.findElementIndex( m_LineA, m_LineType );
  if( iA==-1 )
    {
      Log.LogDebug( "Rule %s: line %s not found.", Name.c_str(), m_LineA.c_str() );
      return;
    }
  if( LineModelElementList[iA]->GetSize()>1 )
    {
      Log.LogDebug( "Rule %s: line %s has size < 1.", Name.c_str(), m_LineA.c_str() );
      iA=-1;
    }
  Int32 iB = LineModelElementList.findElementIndex( m_LineB, m_LineType );
  if( iB==-1 )
    {
      Log.LogDebug( "Rule %s: line %s not found.", Name.c_str(), m_LineB.c_str() );
      return;
    }
  if( LineModelElementList[iB]->GetSize()>1 )
    {
      Log.LogDebug( "Rule %s: line %s has size < 1.", Name.c_str(), m_LineB.c_str() );
      iB=-1;
    }
  if( iA==-1 || iB==-1 || iA==iB )
    {
      Log.LogDebug( "Rule %s: line %s has same index as line %s.", Name.c_str(), m_LineA.c_str(), m_LineB.c_str() );
      return;
    }
  if( LineModelElementList[iA]->IsOutsideLambdaRange() == false )
    {
      Float64 ampA = LineModelElementList[iA]->GetFittedAmplitude( 0 );
      Float64 ampB = LineModelElementList[iB]->GetFittedAmplitude( 0 );

      if( !(ampA<=0.0 && ampB<=0.0) )
      {
          //*
          //Method 0, limit the weakest line's amplitude, no noise taken into account
          Float64 maxB = (m_Coefficient*ampA);
          if(maxB==std::min(maxB, ampB))
          {
              LineModelElementList[iB]->LimitFittedAmplitude(0, maxB);
              //log the correction
              {
                  std::string strTmp0 = boost::str( (boost::format("correct - %-10s") % "2_SINGLE_LINES_AMPLITUDE" ));
                  Logs.append(strTmp0.c_str());
                  std::string strTmp = boost::str( (boost::format("\n\tlineWeak=%-10s, lineStrong=%-10s, previousAmp=%.4e, correctedAmp=%.4e") % m_LineB % m_LineA % ampB % maxB) );
                  Logs.append(strTmp.c_str());
              }
          }
          //*/

          /*
          //Method 1, limit the weakest line's amplitude, only the strongest line's noise is taken into account
          Float64 maxB = (m_Coefficient*ampA) + (erA*nSigma*m_Coefficient);
          if(maxB==std::min(maxB, ampB))
          {
              LineModelElementList[iB]->LimitFittedAmplitude(0, maxB);
              //log the correction
              {
                  std::string strTmp0 = boost::str( (boost::format("correct - %-10s") % "2_SINGLE_LINES_AMPLITUDE" ));
                  Logs.append(strTmp0.c_str());
                  std::string strTmp = boost::str( (boost::format("\n\tlineWeak=%-10s, lineStrong=%-10s, previousAmp=%.4e, correctedAmp=%.4e") % m_LineB % m_LineA % ampB % maxB) );
                  Logs.append(strTmp.c_str());
              }
          }
          //*/

          /*
          //Method 2, correct both lines depending on their sigmas
          if( ampB!=0.0 && (erA!=0 && erB!=0) && std::abs( ampB )>std::abs( ampA*m_Coefficient ) )
          {
              Float64 R = 1.0/m_Coefficient;
              Float64 wA = 0.0;
              if( erA!=0.0 )
              {
                  wA = 1.0/(erA*erA);
              }
              Float64 wB = 0.0;
              if( erB!=0.0 )
              {
                  wB = 1.0/(erB*erB*R*R);
              }
              Float64 correctedA = (ampA*wA + ampB*wB*R)/(wA+wB);
              Float64 correctedB = correctedA/R;
              LineModelElementList[iA]->SetFittedAmplitude( correctedA, erA ); //check: keep the original error sigma ?
              LineModelElementList[iB]->SetFittedAmplitude( correctedB, erB ); //check: keep the original error sigma ?
          }
          else
          {
              if( ampB!=0.0 && ampA==0.0 )
              {
                  Float64 maxB = erA;
                  LineModelElementList[iB]->LimitFittedAmplitude( 0, maxB );
              }
          }
          //*/
      }
  }
}

bool CRule2SingleLinesAmplitude::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

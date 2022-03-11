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
#include "RedshiftLibrary/line/ruleSuperStrongHighest.h"

using namespace NSEpic;
using namespace std;

void CRuleSuperStrong::SetUp( bool EnabledArgument, ... )
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
  va_end(Arguments);
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
  
  for(Int32 iedx = 0; iedx<LineModelElementList.size();iedx++)
  {
    Int32 nLines = LineModelElementList[iedx]->m_Lines.size();
  
    for( Int32 iLineWeak=0; iLineWeak<nLines; iLineWeak++ ) //loop on the weak lines
    {

      bool foundSuperStringTag=false;
      for(Int32 k=0; k<m_SuperStrongTags.size(); k++)
      {
          if(LineModelElementList[iedx]->m_Lines[iLineWeak].GetName() == m_SuperStrongTags[k])
          {
              foundSuperStringTag=true;
              break;
          }
      }
      if( foundSuperStringTag )
      {
          continue;
      }

      if( LineModelElementList[iedx]->m_Lines[iLineWeak].GetType() != m_LineType )
	    {
	      continue;
	    }
      if( LineModelElementList[iedx]->IsOutsideLambdaRange( iLineWeak ) == true )
      {
        continue;
      }
      //Log.LogDebug( "Rule %s: element %d has force weak, type %d and is not outside lambda range.", Name.c_str(), iLineWeak, m_LineType );
      Float64 nSigma = 1.0;
      Float64 ampA = maxiStrong;
      Float64 erA = erStrong;
      Float64 ampB = LineModelElementList[iedx]->GetFittedAmplitude( iLineWeak );

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
          LineModelElementList[iedx]->LimitFittedAmplitude( iLineWeak, maxB );
          //log the correction
          {
              std::string nameWeak = LineModelElementList[iedx]->m_Lines[iLineWeak].GetName();
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
}

bool CRuleSuperStrong::Check( CLineModelElementList& LineModelElementList )
{
  return false;
}

/**
 * \brief Returns the maximum amplitude between strong lines within the support of m_Elements. The referenced er argument will hold the error sigma for the same element.
 **/
Float64 CRuleSuperStrong::FindHighestSuperStrongLineAmp( TStringList superstrongTags , Float64 &er, std::string &name, CLineModelElementList& LineModelElementList )
{
  Float64 maxi = -1.0;
  for(Int32 iedx = 0; iedx<LineModelElementList.size();iedx++)
  {
    Int32 nLines = LineModelElementList[iedx]->m_Lines.size();
    for( Int32 iLineStrong=0; iLineStrong<nLines; iLineStrong++ ) //loop on the strong lines 
    {
        if( LineModelElementList[iedx]->m_Lines[iLineStrong].GetForce() != CLine::nForce_Strong )
        {
          continue;
        }
            if( LineModelElementList[iedx]->m_Lines[iLineStrong].GetType() != m_LineType )
        {
          continue;
        }
        bool foundSuperStringTag=false;
        for(Int32 k=0; k<superstrongTags.size(); k++)
        {
          if(LineModelElementList[iedx]->m_Lines[iLineStrong].GetName() == superstrongTags[k])
          {
              foundSuperStringTag=true;
              break;
          }
        }
        if( !foundSuperStringTag )
        {
          continue;
        }

        if( LineModelElementList[iedx]->IsOutsideLambdaRange( iLineStrong ) == true )
        {
          continue;
        }

        Float64 ampStrong = LineModelElementList[iedx]->GetFittedAmplitude( iLineStrong );
        Float64 erStrong = LineModelElementList[iedx]->GetFittedAmplitudeErrorSigma( iLineStrong );
        // if(erStrong>0.0 && ampStrong>0.0)
        // {
        //     lineSnr = ampStrong/erStrong;
        // }
        if( maxi<ampStrong /*&& lineSnr>validSNRCut*/ )
        {
          maxi = ampStrong;
          er = erStrong;
          name = LineModelElementList[iedx]->m_Lines[iLineStrong].GetName();
        }
      }
  }
  //Log.LogDebug( "Highest strong line amplitude = %f", maxi );
  return maxi;
}

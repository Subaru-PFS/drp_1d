#include <iostream>

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/log/log.h>
#include <RedshiftLibrary/ray/regulament.h>
#include <RedshiftLibrary/ray/linetags.h>

// To be removed once JSON code is in <--
#include <RedshiftLibrary/ray/rule2singlelinesamplitude.h>
#include <RedshiftLibrary/ray/ruleBalmerLinearSolver.h>
#include <RedshiftLibrary/ray/ruleOIIRatioRange.h>
#include <RedshiftLibrary/ray/ruleStrongHigherThanWeak.h>
#include <RedshiftLibrary/ray/ruleSuperStrongHighest.h>
// -->

using namespace NSEpic;
using namespace std;

void CRegulament::Apply( CLineModelElementList& LineModelElementList )
{
    
  for( std::vector<CRule*>::iterator it = m_RulesVector.begin(); it != m_RulesVector.end(); it++ )
    {
      (*it)->Apply ( LineModelElementList );
      std::string logRule = (*it)->GetLogs ( );
      if(logRule.size()>0 && m_LogsEnabled){
        m_RulesLog.push_back(logRule);
      }
  }
}
/*
void CRegulament::ApplyWithRedshift( Float64 Redshift )
{

}
*/
Bool CRegulament::CreateRulesFromJSONFiles( void )
{
  // To be removed once JSON code is in <--
  Bool True = true;
  linetags ltags;
  //CRuleBalmerLinearSolver* ARule8 = new CRuleBalmerLinearSolver( );
  //ARule8->SetUp( True );
  //m_RulesVector.push_back( dynamic_cast<CRule*>( ARule8 ) );

  CRule2SingleLinesAmplitude* ARule1 = new CRule2SingleLinesAmplitude( );
  ARule1->SetUp( True, CRay::nType_Emission, ltags.halpha_em, ltags.hbeta_em, 1.0/2.86*1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule1 ) );
  CRule2SingleLinesAmplitude* ARule2 = new CRule2SingleLinesAmplitude( );
  ARule2->SetUp( True, CRay::nType_Emission, ltags.hbeta_em, ltags.hgamma_em, 0.47*1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule2 ) );
  CRule2SingleLinesAmplitude* ARule3 = new CRule2SingleLinesAmplitude( );
  ARule3->SetUp( True, CRay::nType_Emission, ltags.hgamma_em, ltags.hdelta_em, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule3 ) );
  CRule2SingleLinesAmplitude* ARule4 = new CRule2SingleLinesAmplitude( );
  ARule4->SetUp( True, CRay::nType_Emission, ltags.hdelta_em, ltags.h8_em, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule4 ) );
  CRule2SingleLinesAmplitude* ARule5 = new CRule2SingleLinesAmplitude( );
  ARule5->SetUp( True, CRay::nType_Emission, ltags.h8_em, ltags.h9_em, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule5 ) );
  CRule2SingleLinesAmplitude* ARule6 = new CRule2SingleLinesAmplitude( );
  ARule6->SetUp( True, CRay::nType_Emission, ltags.h9_em, ltags.h10_em, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule6 ) );
  CRule2SingleLinesAmplitude* ARule7 = new CRule2SingleLinesAmplitude( );
  ARule7->SetUp( True, CRay::nType_Emission, ltags.h10_em, ltags.h11_em, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule7 ) );

  CRule2SingleLinesAmplitude* ARule2A = new CRule2SingleLinesAmplitude( );
  ARule2A->SetUp( True, CRay::nType_Emission, ltags.hbeta_abs, ltags.hgamma_abs, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule2A ) );
  CRule2SingleLinesAmplitude* ARule3A = new CRule2SingleLinesAmplitude( );
  ARule3A->SetUp( True, CRay::nType_Emission, ltags.hgamma_abs, ltags.hdelta_abs, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule3A ) );
  CRule2SingleLinesAmplitude* ARule4A = new CRule2SingleLinesAmplitude( );
  ARule4A->SetUp( True, CRay::nType_Emission, ltags.hdelta_abs, ltags.h8_abs, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule4A ) );
  CRule2SingleLinesAmplitude* ARule5A = new CRule2SingleLinesAmplitude( );
  ARule5A->SetUp( True, CRay::nType_Emission, ltags.h8_abs, ltags.h9_abs, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule5A ) );
  CRule2SingleLinesAmplitude* ARule6A = new CRule2SingleLinesAmplitude( );
  ARule6A->SetUp( True, CRay::nType_Emission, ltags.h9_abs, ltags.h10_abs, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule6A ) );
  CRule2SingleLinesAmplitude* ARule7A = new CRule2SingleLinesAmplitude( );
  ARule7A->SetUp( True, CRay::nType_Emission, ltags.h10_abs, ltags.h11_abs, 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule7A ) );


  CRuleRatioRange* ARule9 = new CRuleRatioRange( );
  ARule9->SetUp( True, CRay::nType_Emission, ltags.oII3726_em, ltags.oII3729_em, 2.5 );

  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule9 ) );
  CRuleStrongHigherThanWeak* ARule10 = new CRuleStrongHigherThanWeak( );
  ARule10->SetUp( True, CRay::nType_Emission );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule10 ) );
  CRuleStrongHigherThanWeak* ARule11 = new CRuleStrongHigherThanWeak( );
  ARule11->SetUp( True, CRay::nType_Absorption );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule11 ) );


  CRuleRatioRange* ARule13 = new CRuleRatioRange( );
  ARule13->SetUp( True, CRay::nType_Emission, ltags.cIII1907_em, ltags.cIII1909_em, 2.0 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule13 ) );

  // OII and Halpha Super Strong
  //CRuleSuperStrong* ARule12 = new CRuleSuperStrong( );
  //ARule12->SetUp( True, CRay::nType_Emission, ltags.oII3726_em, ltags.oII3729_em, ltags.halpha_em, 1.1 );
  //m_RulesVector.push_back( dynamic_cast<CRule*>( ARule12 ) );

  if(m_LogsEnabled)
  {
    m_RulesLog.push_back("Linemodel-Regulament: rules creation");
  }
  return true;
}

void CRegulament::EnableRulesAccordingToParameters( std::string Parameters )
{
  if( Parameters=="no" )
    {
        return;
    }
  for( std::vector<CRule*>::iterator it = m_RulesVector.begin(); it != m_RulesVector.end(); it++ )
    {
      Bool enableRule = Parameters.find ( (*it)->Name ) != std::string::npos;
      if( Parameters=="all" || enableRule )
    {
      Log.LogDebug( "Enabling rule %s.", (*it)->Name.c_str() );
      (*it)->Enabled = true;
    }else
      {
          Log.LogDebug( "Disabling rule %s.", (*it)->Name.c_str() );
          (*it)->Enabled = false;
      }
    }
}

void CRegulament::EnableLogs( bool enable )
{
    m_LogsEnabled = enable;
    m_RulesLog.clear();
}
std::vector<string> CRegulament::GetLogs( )
{
    return m_RulesLog;
}


#include <iostream>

#include <epic/core/common/datatypes.h>
#include <epic/core/log/log.h>
#include <epic/redshift/ray/regulament.h>
// To be removed once JSON code is in <--
#include <epic/redshift/ray/rule2singlelinesamplitude.h>
#include <epic/redshift/ray/ruleBalmerLinearSolver.h>
#include <epic/redshift/ray/ruleOIIRatioRange.h>
#include <epic/redshift/ray/ruleStrongHigherThanWeak.h>
#include <epic/redshift/ray/ruleSuperStrongHighest.h>
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

  CRuleBalmerLinearSolver* ARule8 = new CRuleBalmerLinearSolver( );
  ARule8->SetUp( True );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule8 ) );

  CRule2SingleLinesAmplitude* ARule1 = new CRule2SingleLinesAmplitude( );
  ARule1->SetUp( True, CRay::nType_Emission, std::string( "Halpha" ).c_str(), std::string( "Hbeta" ).c_str(), 1.0/2.86*1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule1 ) );
  CRule2SingleLinesAmplitude* ARule2 = new CRule2SingleLinesAmplitude( );
  ARule2->SetUp( True, CRay::nType_Emission, std::string( "Hbeta" ).c_str(), std::string( "Hgamma" ).c_str(), 0.47*1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule2 ) );
  CRule2SingleLinesAmplitude* ARule3 = new CRule2SingleLinesAmplitude( );
  ARule3->SetUp( True, CRay::nType_Emission, std::string( "Hgamma" ).c_str(), std::string( "Hdelta" ).c_str(), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule3 ) );
  CRule2SingleLinesAmplitude* ARule4 = new CRule2SingleLinesAmplitude( );
  ARule4->SetUp( True, CRay::nType_Emission, std::string( "Hdelta" ).c_str(), std::string( "H8" ).c_str(), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule4 ) );
  CRule2SingleLinesAmplitude* ARule5 = new CRule2SingleLinesAmplitude( );
  ARule5->SetUp( True, CRay::nType_Emission, std::string( "H8" ).c_str(), std::string( "H9" ).c_str(), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule5 ) );
  CRule2SingleLinesAmplitude* ARule6 = new CRule2SingleLinesAmplitude( );
  ARule6->SetUp( True, CRay::nType_Emission, std::string( "H9" ).c_str(), std::string( "H10" ).c_str(), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule6 ) );
  CRule2SingleLinesAmplitude* ARule7 = new CRule2SingleLinesAmplitude( );
  ARule7->SetUp( True, CRay::nType_Emission, std::string( "H10" ).c_str(), std::string( "H11" ).c_str(), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule7 ) );

  CRuleOIIRatioRange* ARule9 = new CRuleOIIRatioRange( );
  ARule9->SetUp( True, CRay::nType_Emission, "[OII]3726", "[OII]3729", 2.0 );

  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule9 ) );
  CRuleStrongHigherThanWeak* ARule10 = new CRuleStrongHigherThanWeak( );
  ARule10->SetUp( True, CRay::nType_Emission );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule10 ) );
  CRuleStrongHigherThanWeak* ARule11 = new CRuleStrongHigherThanWeak( );
  ARule11->SetUp( True, CRay::nType_Absorption );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule11 ) );


  CRuleOIIRatioRange* ARule13 = new CRuleOIIRatioRange( );
  ARule13->SetUp( True, CRay::nType_Emission, "[CIII]1907", "[CIII]1909", 2.0 );
  m_RulesVector.push_back( dynamic_cast<CRule*>( ARule13 ) );

  // OII and Halpha Super Strong
  //CRuleSuperStrong* ARule12 = new CRuleSuperStrong( );
  //ARule12->SetUp( True, CRay::nType_Emission, std::string( "[OII]3726" ).c_str(), std::string( "[OII]3726" ).c_str(), std::string( "Halpha" ).c_str(), 1.1 );
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


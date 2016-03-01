#include <iostream>

#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/regulament.h>
// To be removed once JSON code is in <--
#include <epic/redshift/ray/rule2singlelinesamplitude.h>
// -->

using namespace NSEpic;
using namespace std;

void CRegulament::Apply( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements )
{
  m_Elements = LinemodelElements;
  for( std::vector<CRule*>::iterator it = m_RulesVector.begin(); it != m_RulesVector.end(); it++ )
    {
      (*it)->Apply ( LinemodelElements );
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
  CRule2SingleLinesAmplitude* ARule1 = new CRule2SingleLinesAmplitude( );
  ARule1->SetUp ( True, CRay::nType_Emission, std::string( "Halpha" ), std::string( "Hbeta" ), 1.0/2.86*1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*> ( ARule1 ) );
  CRule2SingleLinesAmplitude* ARule2 = new CRule2SingleLinesAmplitude( );
  ARule2->SetUp ( True, CRay::nType_Emission, std::string( "Hbeta" ), std::string( "Hgamma" ), 0.47*1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*> ( ARule2 ) );
  CRule2SingleLinesAmplitude* ARule3 = new CRule2SingleLinesAmplitude( );
  ARule3->SetUp ( True, CRay::nType_Emission, std::string( "Hgamma" ), std::string( "Hdelta" ), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*> ( ARule3 ) );
  CRule2SingleLinesAmplitude* ARule4 = new CRule2SingleLinesAmplitude( );
  ARule4->SetUp ( True, CRay::nType_Emission, std::string( "Hdelta" ), std::string( "H8" ), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*> ( ARule4 ) );
  CRule2SingleLinesAmplitude* ARule5 = new CRule2SingleLinesAmplitude( );
  ARule5->SetUp ( True, CRay::nType_Emission, std::string( "H8" ), std::string( "H9" ), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*> ( ARule5 ) );
  CRule2SingleLinesAmplitude* ARule6 = new CRule2SingleLinesAmplitude( );
  ARule6->SetUp ( True, CRay::nType_Emission, std::string( "H9" ), std::string( "H10" ), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*> ( ARule6 ) );
  CRule2SingleLinesAmplitude* ARule7 = new CRule2SingleLinesAmplitude( );
  ARule7->SetUp ( True, CRay::nType_Emission, std::string( "H10" ), std::string( "H11" ), 1.1 );
  m_RulesVector.push_back( dynamic_cast<CRule*> ( ARule7 ) );
  // -->
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
	  (*it)->Enabled = true;
	}
    }
}

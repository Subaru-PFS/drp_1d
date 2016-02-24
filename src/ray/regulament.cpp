#include <epic/redshift/ray/regulament.h>

using namespace NSEpic;
using namespace std;

void CRegulament::Apply( std::vector<boost::shared_ptr<CLineModelElement> >& LinemodelElements )
{
  m_Elements = LinemodelElements;
  for( std::vector<CRule>::iterator it = m_RulesVector.begin(); it != m_RulesVector.end(); it++ )
    {
      (*it).Apply ( LinemodelElements );
    }
}
/*
void CRegulament::ApplyWithRedshift( Float64 Redshift )
{
  
}
*/
Bool CRegulament::CreateRulesFromJSONFiles( void )
{
  return false;
}

void CRegulament::EnableRulesAccordingToParameters( std::string Parameters )
{
  if( Parameters=="no" )
    {
        return;
    }
  for( std::vector<CRule>::iterator it = m_RulesVector.begin(); it != m_RulesVector.end(); it++ )
    {
      Bool enableRule = Parameters.find ( (*it).Name ) != std::string::npos;
      if( Parameters=="all" || enableRule )
	{
	  (*it).Enabled = true;
	}
    }
}

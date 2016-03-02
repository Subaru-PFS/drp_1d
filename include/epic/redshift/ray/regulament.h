#ifndef _REDSHIFT_RAY_REGULAMENT_
#define _REDSHIFT_RAY_REGULAMENT_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/linemodel/elementlist.h>
#include <epic/redshift/ray/rule.h>

#include <boost/format.hpp>

namespace NSEpic
{
  /**
   * \ingroup Redshift
   * \brief Control class for preparing and applying Linemodel rules.
   */
  class CRegulament
  {
  public:
    void Apply( CLineModelElementList& LineModelElementList );
    //void ApplyWithRedshift( std::vector<boost::shared_ptr<CLineModelElement> > LinemodelElements, Float64 Redshift );
    Bool CreateRulesFromJSONFiles( void );
    void EnableRulesAccordingToParameters( std::string Parameters );
  private:
    Float64 m_Redshift;
    std::vector<CRule*> m_RulesVector;
  };
}

#endif

#ifndef _REDSHIFT_RAY_RULE_
#define _REDSHIFT_RAY_RULE_

#include <epic/core/common/datatypes.h>
#include <boost/format.hpp>
#include <epic/redshift/linemodel/elementlist.h>

namespace NSEpic
{
  /**
   * \ingroup Redshift
   * Abstract class for common functionality of rules.
   */
  class CRule
  {
  public:
    Bool Enabled;
    CRule ( );
    void Apply( std::vector<boost::shared_ptr<CLineModelElement> > LinemodelElements );
    virtual Bool Check( std::vector<boost::shared_ptr<CLineModelElement> > LinemodelElements ) = 0;
  private:
    virtual void Correct( std::vector<boost::shared_ptr<CLineModelElement> > LinemodelElements ) = 0;
  };
}

#endif


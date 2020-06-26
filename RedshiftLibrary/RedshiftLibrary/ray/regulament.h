#ifndef _REDSHIFT_RAY_REGULAMENT_
#define _REDSHIFT_RAY_REGULAMENT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/linemodel/elementlist.h>
#include <RedshiftLibrary/ray/rule.h>

#include <boost/format.hpp>

namespace NSEpic {
/**
 * \ingroup Redshift
 * \brief Control class for preparing and applying Linemodel rules.
 */
class CRegulament
{
  public:
    CRegulament();
    ~CRegulament();

    void Apply(CLineModelElementList &LineModelElementList);
    // void ApplyWithRedshift( std::vector<boost::shared_ptr<CLineModelElement>
    // > LinemodelElements, Float64 Redshift );
    Bool CreateRulesFromJSONFiles(void);
    void EnableRulesAccordingToParameters(std::string Parameters);
    std::vector<std::string> GetLogs();
    void EnableLogs(bool enable);

  private:
    Float64 m_Redshift;
    std::vector<CRule *> m_RulesVector;
    std::vector<std::string> m_RulesLog;
    bool m_LogsEnabled = false;
};
} // namespace NSEpic

#endif

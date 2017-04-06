#ifndef _REDSHIFT_RAY_CATALOGSOFFSETS_
#define _REDSHIFT_RAY_CATALOGSOFFSETS_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/linemodel/elementlist.h>


#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * /ingroup Redshift

 */
class CRayCatalogsOffsets
{

public:
    CRayCatalogsOffsets();
    ~CRayCatalogsOffsets();
    Bool Init(std::string calibrationPath);
    Bool SetCtlgRelPath( const char* relPath );

    Bool Load( const char* filePath );
    Bool SetLinesOffsets(CLineModelElementList &LineModelElementList);

private:
    std::string m_Catalog_relpath;
    std::vector<Float64> m_Offsets;
    std::vector<std::string> m_Names;

};


}

#endif

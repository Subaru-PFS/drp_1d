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
class CLineCatalogsOffsets
{

public:

    struct SOffsetsCatalog
    {
        std::string filePath;
        std::vector<Float64> Offsets;
        std::vector<std::string> Names;
    };

    CLineCatalogsOffsets();
    ~CLineCatalogsOffsets();
    Bool Init(std::string calibrationPath);
    Bool SetCtlgRelPath( const char* relPath );

    Bool Load( const char* dirPath );
    Bool LoadCatalog( const char* filePath );
    Bool SetLinesOffsets(CLineModelElementList &LineModelElementList, Int32 index);

    // Hack: select/bypass stack automatically from its name in the reference_stack catalog
    Bool SetLinesOffsetsAutoSelectStack(CLineModelElementList &LineModelElementList, std::string spectrumName);
    Int32 AutoSelectStackFromReferenceFile(std::string spectrumName);

private:
    std::string m_Catalogs_relpath;
    std::string m_Calibration_path;

    std::vector<SOffsetsCatalog> m_OffsetsCatalog;

};


}

#endif

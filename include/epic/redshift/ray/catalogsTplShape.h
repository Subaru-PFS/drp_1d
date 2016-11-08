#ifndef _REDSHIFT_RAY_CATALOGSTPLSHAPE_
#define _REDSHIFT_RAY_CATALOGSTPLSHAPE_

#include <epic/core/common/datatypes.h>
#include <epic/redshift/ray/ray.h>
#include <epic/redshift/ray/catalog.h>


#include <boost/format.hpp>

#include <vector>
#include <string>

namespace NSEpic
{

/**
 * /ingroup Redshift

 */
class CRayCatalogsTplShape
{

public:
    CRayCatalogsTplShape();
    ~CRayCatalogsTplShape();
    Bool Init();
    Bool Load( const char* dirPath );
    //Bool AreCatalogsAligned( const CRayCatalog::TRayVector& restRayList, Int32 typeFilter, Int32 forceFilter  );
    Float64 GetBestFit(const CRayCatalog::TRayVector& restRayList, std::vector<Float64> fittedAmplitudes, std::vector<Float64> fittedErrors, std::vector<Float64> &amplitudesCorrected , std::__cxx11::string &bestTplName);

private:
    Float64 GetFit(std::vector<Float64> ampsLM, std::vector<Float64> errLM, std::vector<Float64> ampsTPL , std::vector<Float64> &ampsCorrected);


    std::vector<std::string> m_RayCatalogNames;
    std::vector<CRayCatalog> m_RayCatalogList;

};


}

#endif

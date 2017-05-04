#ifndef _REDSHIFT_RAY_CATALOGSTPLSHAPE_
#define _REDSHIFT_RAY_CATALOGSTPLSHAPE_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/ray.h>
#include <RedshiftLibrary/ray/catalog.h>
#include <RedshiftLibrary/linemodel/elementlist.h>


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
    Bool Init(std::string calibrationPath);
    Bool SetTplctlgRelPath( const char* relPath );

    Bool Load( const char* dirPath );
    bool LoadVelocities( const char* filepath, Int32 k );
    //Bool AreCatalogsAligned( const CRayCatalog::TRayVector& restRayList, Int32 typeFilter, Int32 forceFilter  );
    Float64 GetBestFit(const CRayCatalog::TRayVector& restRayList, std::vector<Float64> fittedAmplitudes, std::vector<Float64> fittedErrors, std::vector<Float64> &amplitudesCorrected , std::string &bestTplName);
    CRayCatalog::TRayVector GetRestLinesList( const Int32 index );
    Int32 GetCatalogsCount();
    std::string GetCatalogName(Int32 idx);
    Bool GetCatalogVelocities(Int32 idx, Float64& elv, Float64& alv );
    Bool SetMultilineNominalAmplitudes(CLineModelElementList& LineModelElementList, Int32 iLine);
    Bool SetLyaProfile(CLineModelElementList &LineModelElementList, Int32 iCatalog);
    Bool InitLineCorrespondingAmplitudes(CLineModelElementList &LineModelElementList);
    Bool SetMultilineNominalAmplitudesFast(CLineModelElementList &LineModelElementList, Int32 iCatalog);

private:
    Float64 GetFit(std::vector<Float64> ampsLM, std::vector<Float64> errLM, std::vector<Float64> ampsTPL , std::vector<Float64> &ampsCorrected);

    std::string tplshapedcatalog_relpath;

    std::vector<std::string> m_RayCatalogNames;
    std::vector<CRayCatalog> m_RayCatalogList;
    std::vector<std::vector<Float64>> m_RayCatalogLinesCorrespondingNominalAmp;
    std::vector<Float64> m_ELvelocities;
    std::vector<Float64> m_ABSvelocities;

};


}

#endif

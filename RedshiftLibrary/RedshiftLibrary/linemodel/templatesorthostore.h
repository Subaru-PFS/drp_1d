#ifndef _REDSHIFT_LINEMODEL_TEMPLATES_ORTHO_STORE_
#define _REDSHIFT_LINEMODEL_TEMPLATES_ORTHO_STORE_


#include <epic/core/common/datatypes.h>
#include <epic/redshift/spectrum/template/template.h>
#include <epic/redshift/spectrum/template/catalog.h>


#include <boost/filesystem.hpp>
#include <vector>

namespace NSEpic
{

class CTemplatesOrthoStore
{
public:

    struct SOrthoFitting{
        std::vector<Float64>    lambda;
        std::vector<Float64>    mtmCumulative;
    };

    struct SCatalogDescription{
        Float64                 velocityEmission;
        Float64                 velocityAbsorption;
        //std::string             profile;
    };


    typedef std::vector< CTemplatesOrthoStore::SOrthoFitting >          TTemplatesOrthoFittingRefList;
    typedef std::map< std::string, TTemplatesOrthoFittingRefList >          TTemplatesOrthoFittingRefDict;


    CTemplatesOrthoStore();
    ~CTemplatesOrthoStore();
    bool Add(std::shared_ptr<CTemplateCatalog> tplCtlg);
    std::shared_ptr<CTemplateCatalog> getTplCatalog(Int32 ctlgIdx);

private:
    std::vector<std::shared_ptr<CTemplateCatalog>>   m_CatalogList;
    std::vector<TTemplatesOrthoFittingRefDict>    m_OrthoFittingList;
    std::vector<SCatalogDescription>    m_DescriptionList;

};



}

#endif

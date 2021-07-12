#include "RedshiftLibrary/linemodel/templatesorthostore.h"
#include "RedshiftLibrary/linemodel/elementlist.h"



using namespace NSEpic;


CTemplatesOrthoStore::CTemplatesOrthoStore()
{

}

/**
 * \brief Empty destructor.
 **/
CTemplatesOrthoStore::~CTemplatesOrthoStore()
{

}


bool CTemplatesOrthoStore::Add(std::shared_ptr<CTemplateCatalog> tplCtlg)
{
    m_CatalogList.push_back(tplCtlg);
    return 0;
}


std::shared_ptr<CTemplateCatalog> CTemplatesOrthoStore::getTplCatalog(Int32 ctlgIdx)
{
    if(ctlgIdx>=m_CatalogList.size())
    {
        return NULL;
    }
    return m_CatalogList[ctlgIdx];
}

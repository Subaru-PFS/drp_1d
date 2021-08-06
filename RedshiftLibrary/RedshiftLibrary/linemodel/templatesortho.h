#ifndef _REDSHIFT_LINEMODEL_TEMPLATESORTHO_
#define _REDSHIFT_LINEMODEL_TEMPLATESORTHO_

#include "RedshiftLibrary/common/range.h"
#include "RedshiftLibrary/common/datatypes.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "RedshiftLibrary/ray/catalog.h"
#include "RedshiftLibrary/spectrum/spectrum.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/operator/templatefitting.h"

#include "RedshiftLibrary/operator/linemodelresult.h"
#include "RedshiftLibrary/linemodel/element.h"

#include "RedshiftLibrary/spectrum/template/catalog.h"
#include "RedshiftLibrary/linemodel/templatesorthostore.h"

#include <boost/shared_ptr.hpp>

#include <memory>


namespace NSEpic
{  
class CTemplatesOrthogonalization
{

public:

    CTemplatesOrthogonalization(const CTemplateCatalog &tplCatalog,
                                const TStringList &tplCategoryList,
                                const std::string calibrationPath,
                                const CRayCatalog::TRayVector &restRayList,
                                const std::string &opt_fittingmethod,
                                const std::string &widthType,
                                const Float64 opt_nsigmasupport,
                                const Float64 velocityEmission,
                                const Float64 velocityAbsorption,
                                const std::string &opt_rules,
                                const std::string &opt_rigidity,
                                std::shared_ptr<const CLSF> lsf,
                                bool enableOrtho=false);

    ~CTemplatesOrthogonalization();

    CTemplateCatalog getOrthogonalTplCatalog();
    CTemplatesOrthoStore getOrthogonalTplStore();

private:

    bool m_enableOrtho;
    std::shared_ptr<const CLSF> m_LSF = nullptr;
    Int32 OrthogonalizeTemplate(const CTemplate& inputTemplate,
                                const std::string calibrationPath,
                                const CRayCatalog::TRayVector &restRayList,
                                const std::string &opt_fittingmethod,
                                const std::string &widthType,
                                const Float64 opt_nsigmasupport,
                                const Float64 velocityEmission,
                                const Float64 velocityAbsorption,
                                const std::string &opt_rules,
                                const std::string &opt_rigidity);


    CTemplateCatalog m_tplCatalogOrthogonal = CTemplateCatalog("zero"); //note: no need to estimate continuum free templates here, //note2: todo: bound to disappear when the tplorthostore is fully implemented
    CTemplatesOrthoStore m_tplOrthoStore;

};

}

#endif // _REDSHIFT_LINEMODEL_TEMPLATESORTHO_

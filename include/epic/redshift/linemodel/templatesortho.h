#ifndef TEMPLATESORTHO_H
#define TEMPLATESORTHO_H

#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <epic/redshift/operator/chisquare2.h>

#include <epic/redshift/operator/linemodelresult.h>
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/singleline.h>

#include <epic/redshift/spectrum/template/catalog.h>

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
                                const std::string &opt_continuumcomponent,
                                const std::string &widthType,
                                const Float64 resolution,
                                const Float64 velocityEmission,
                                const Float64 velocityAbsorption,
                                const std::string &opt_rules,
                                const std::string &opt_rigidity);

    ~CTemplatesOrthogonalization();

    CTemplateCatalog getOrthogonalTplCatalog();

private:

    Int32 OrthogonalizeTemplate(const CTemplate& inputTemplate,
                                const std::string calibrationPath,
                                const CRayCatalog::TRayVector &restRayList,
                                const std::string &opt_fittingmethod,
                                const std::string &widthType,
                                const Float64 resolution,
                                const Float64 velocityEmission,
                                const Float64 velocityAbsorption,
                                const std::string &opt_rules,
                                const std::string &opt_rigidity);


    CTemplateCatalog m_tplCatalogOrthogonal;

};

}







#endif // TEMPLATESORTHO_H

#ifndef ELEMENTLIST_H
#define ELEMENTLIST_H


#include <epic/core/common/ref.h>
#include <epic/core/common/range.h>
#include <epic/redshift/common/datatypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <epic/redshift/ray/catalog.h>
#include <epic/redshift/spectrum/spectrum.h>

#include <epic/redshift/operator/linemodelresult.h>
#include <epic/redshift/linemodel/element.h>
#include <epic/redshift/linemodel/singleline.h>
#include <boost/shared_ptr.hpp>

namespace NSEpic
{

class CLineModelElementList
{

public:

    CLineModelElementList( const CSpectrum& spectrum, const CRayCatalog::TRayVector& restRayList );
    ~CLineModelElementList();

    void LoadCatalog(const CRayCatalog::TRayVector& restRayList);
    void LoadCatalogMultilineBalmer(const CRayCatalog::TRayVector& restRayList);
    void LoadCatalogSingleLines(const CRayCatalog::TRayVector& restRayList);

    void fit(Float64 redshift, CLineModelResult::SLineModelSolution &modelSolution);
    void addToModel();

    CLineModelResult::SLineModelSolution GetModelSolution();
    const CSpectrum&                GetModelSpectrum() const;

private:

    std::vector<Int32> findLineIdxInCatalog(const CRayCatalog::TRayVector& restRayList, std::string strTag);
    Void Apply2SingleLinesAmplitudeRule(std::string lineA, std::string lineB, Float64 coeff );

    void addSingleLine(const CRay &r, Int32 index, Float64 nominalWidth);
    void addDoubleLine(const CRay &r1, const CRay &r2, Int32 index1, Int32 index2, Float64 nominalWidth, Float64 a1, Float64 a2);

    void applyRules();
    Int32 FindElementIndex(Int32 LineCatalogIndex);
    Int32 FindElementIndex(std::string LineTagStr);


    std::vector<boost::shared_ptr<CLineModelElement>  > m_Elements;
    //std::vector<CLineModelElement*> m_Elements;
    //std::vector<CSingleLine*> m_Elements;

    CRef<CSpectrum>   m_SpectrumModel;
    CSpectrumFluxAxis m_SpcFluxAxis;

    CRayCatalog::TRayVector m_RestRayList;
};

}







#endif // ELEMENTLIST_H


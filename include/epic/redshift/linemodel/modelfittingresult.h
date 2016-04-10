#ifndef _REDSHIFT_LINEMODEL_MODELFITTINGRESULT_
#define _REDSHIFT_LINEMODEL_MODELFITTINGRESULT_


#include <epic/redshift/processflow/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

#include <epic/redshift/spectrum/spectrum.h>
#include <epic/redshift/operator/linemodelresult.h>


namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
class CModelFittingResult : public COperatorResult
{

public:

    CModelFittingResult(CLineModelResult::SLineModelSolution _lineModelSolution, Float64 _redshift, Float64 _merit, CRayCatalog::TRayVector _restRayList, Float64 _velEmission=-1.0, Float64 _velAbsorption=-1.0 );
    CModelFittingResult();
    virtual ~CModelFittingResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Void Load( const char* filePath );

    const CLineModelResult::SLineModelSolution& GetLineModelSolution() const;


private:

    CLineModelResult::SLineModelSolution LineModelSolution;
    Float64 Redshift;
    Float64 Merit;

    CRayCatalog::TRayVector restRayList;
    Float64 VelocityEmission;
    Float64 VelocityAbsorption;
};

inline
const CLineModelResult::SLineModelSolution& CModelFittingResult::GetLineModelSolution() const
{
    return LineModelSolution;
}


}

#endif

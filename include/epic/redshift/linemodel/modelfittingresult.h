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

    CModelFittingResult( CLineModelResult::SLineModelSolution _lineModelSolution, Float64 _redshift, Float64 _merit, CRayCatalog::TRayVector _restRayList );
    CModelFittingResult();
    virtual ~CModelFittingResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;

private:

    CLineModelResult::SLineModelSolution LineModelSolution;
    Float64 Redshift;
    Float64 Merit;

    CRayCatalog::TRayVector restRayList;

};


}

#endif

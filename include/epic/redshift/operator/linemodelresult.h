#ifndef _REDSHIFT_OPERATOR_LINEMODELRESULT_
#define _REDSHIFT_OPERATOR_LINEMODELRESULT_

#include <epic/redshift/operator/result.h>
#include <epic/core/common/datatypes.h>
#include <epic/redshift/operator/operator.h>

#include <epic/redshift/ray/catalog.h>

namespace NSEpic
{

class CLineModelResult : public COperatorResult
{

    DEFINE_MANAGED_OBJECT( CLineModelResult )

public:
    struct SLineModelSolution
    {
        std::vector<Float64> Amplitudes;
        std::vector<Float64> Widths;
        std::vector<Bool> OutsideLambdaRange;
        std::vector<TInt32Range> fittingIndexRange;
    };

    CLineModelResult();
    virtual ~CLineModelResult();

    Void Save( const COperatorResultStore& store, std::ostream& stream ) const;
    Void SaveLine( const COperatorResultStore& store, std::ostream& stream ) const;
    Void Load( std::istream& stream );

    TFloat64List            Redshifts;
    TFloat64List            ChiSquare;
    TFloat64List            Area;
    TFloat64List            Extrema;
    COperator::TStatusList  Status;

    std::vector<SLineModelSolution> LineModelSolutions;
    CRayCatalog::TRayVector restRayList;


};


}

#endif

#ifndef _REDSHIFT_LINEMODEL_MODELCONTINUUMFITTINGRESULT_
#define _REDSHIFT_LINEMODEL_MODELCONTINUUMFITTINGRESULT_


#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/operator/operator.h>

#include <RedshiftLibrary/spectrum/spectrum.h>
#include <RedshiftLibrary/operator/linemodelresult.h>


namespace NSEpic
{

  /**
   * \ingroup Redshift
   */
class CModelContinuumFittingResult : public COperatorResult
{

public:

    CModelContinuumFittingResult(Float64 _redshift,
                                 std::string _name,
                                 Float64 _merit,
                                 Float64 _amp,
                                 Float64 _ismCoeff,
                                 Int32 _igmIndex,
                                 Float64 _fitting_snr);
    CModelContinuumFittingResult();
    virtual ~CModelContinuumFittingResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }


private:

    Float64 Redshift;
    std::string Name;
    Float64 Merit;
    Float64 Amp;
    Float64 IsmCoeff;
    Int32   IgmIndex;

    //fitting info
    Float64 Fitting_snr;
};

}

#endif


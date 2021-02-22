#ifndef _REDSHIFT_LINEMODEL_MODELFITTINGRESULT_
#define _REDSHIFT_LINEMODEL_MODELFITTINGRESULT_


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
class CModelFittingResult : public COperatorResult
{

public:

    CModelFittingResult(CLineModelSolution _lineModelSolution, Float64 _redshift, Float64 _merit, CRayCatalog::TRayVector _restRayList, Float64 _velEmission=-1.0, Float64 _velAbsorption=-1.0 );
    CModelFittingResult();
    virtual ~CModelFittingResult();

    void Save( const CDataStore& store, std::ostream& stream ) const;
    void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    inline Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const
    {
        return 1;
    }
    void Load( const char* filePath );

    const CLineModelSolution& GetLineModelSolution() const;
  void getData(const std::string& name, double **data, int *size) const;
  void getData(const std::string& name, std::string *data, int *size) const;
  void getData(const std::string& name, int  **data, int *size) const;


private:

    CLineModelSolution LineModelSolution;
    Float64 Redshift;
    Float64 Merit;

    CRayCatalog::TRayVector restRayList;
    Float64 VelocityEmission;
    Float64 VelocityAbsorption;

  //copies from restRayList For output only
    mutable std::vector<Int32> rayId;
    mutable std::vector<Float64> rayLambdaRest;

};

inline
const CLineModelSolution& CModelFittingResult::GetLineModelSolution() const
{
    return LineModelSolution;
}


}

#endif

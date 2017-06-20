#ifndef _REDSHIFT_OPERATOR_LINEMODELSOLVERESULT_
#define _REDSHIFT_OPERATOR_LINEMODELSOLVERESULT_

#include <RedshiftLibrary/processflow/result.h>
#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/ray/catalog.h>

#include <vector>

namespace NSEpic
{

class CProcessFlowContext;

/**
 * \ingroup Redshift
 */
class CLineModelSolveResult : public COperatorResult
{

public:

    enum EType
    {
         nType_raw = 1,
         nType_continuumOnly = 2,
         nType_noContinuum = 3,
         nType_all = 4,
    };

    CLineModelSolveResult();
    virtual ~CLineModelSolveResult();

    Void Save( const CDataStore& store, std::ostream& stream ) const;
    Void SaveLine( const CDataStore& store, std::ostream& stream ) const;
    Bool GetBestRedshift(const CDataStore& store, Float64& redshift, Float64& merit , Float64 &sigma) const;
    Bool GetBestRedshiftLogArea( const CDataStore& store, Float64& redshift, Float64& merit ) const;
    Bool GetBestRedshiftWithStrongELSnrPrior( const CDataStore& store, Float64& redshift, Float64& merit ) const;


    Void SetReliabilityLabel( std::string lbl );


private:
    std::string m_ReliabilityLabel;


};


}

#endif


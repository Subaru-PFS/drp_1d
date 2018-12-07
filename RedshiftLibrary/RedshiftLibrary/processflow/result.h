#ifndef _REDSHIFT_OPERATOR_RESULT_
#define _REDSHIFT_OPERATOR_RESULT_

#include <RedshiftLibrary/common/datatypes.h>

#include <vector>
#include <ostream>
#include <map>

namespace NSEpic
{

class CDataStore;

/**
 * \ingroup Redshift
 */
class COperatorResult
{

public:

    COperatorResult();
    virtual ~COperatorResult();

    virtual void Save( const CDataStore& store, std::ostream& stream ) const = 0;
    virtual void SaveLine( const CDataStore& store, std::ostream& stream ) const = 0;
    //virtual void Load( std::istream& stream ) = 0;

    void SetReliabilityLabel( std::string lbl );
    void SetTypeLabel( std::string lbl );
    virtual Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const = 0;


protected:

    std::string m_ReliabilityLabel="-1";
    std::string m_TypeLabel="-1";
};

typedef std::vector< std::shared_ptr<COperatorResult> >           TOperatorResultList;
typedef std::map< std::string, std::shared_ptr< const COperatorResult> > TOperatorResultMap;

}

#endif

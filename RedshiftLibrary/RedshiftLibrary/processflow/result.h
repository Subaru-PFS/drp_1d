#ifndef _REDSHIFT_OPERATOR_RESULT_
#define _REDSHIFT_OPERATOR_RESULT_

#include <RedshiftLibrary/common/datatypes.h>
#include <RedshiftLibrary/continuum/indexes.h>
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
    virtual void SaveJSON(const CDataStore& store, std::ostream& stream) const;
    //virtual void Load( std::istream& stream ) = 0;

    void SetReliabilityLabel( std::string lbl );
    void SetTypeLabel( std::string lbl );
    virtual Int32 GetEvidenceFromPdf(const CDataStore& store, Float64 &evidence) const = 0;

    void SaveFloat64(std::ostream& stream,Float64 data) const;
    void SaveTFloat64List(std::ostream& stream,std::string name,TFloat64List data, TFloat64List order) const;
    void SaveTFloat64ListOfList(std::ostream& stream,std::string name,std::vector<TFloat64List> data, TFloat64List order) const;
    void SaveInt32Vector(std::ostream& stream,std::string name,std::vector<Int32> data, TFloat64List order) const;
    void SaveStringVector(std::ostream& stream,std::string name,std::vector<std::string>, TFloat64List order) const;
    void SaveStringVectorOfVector(std::ostream& stream,std::string name,std::vector<std::vector<std::string>>, TFloat64List order) const;
    void SaveTContinuumIndexListVector(std::ostream& stream,std::string name,std::vector<CContinuumIndexes::TContinuumIndexList>, TFloat64List order) const;
protected:

    std::string m_ReliabilityLabel="-1";
    std::string m_TypeLabel="-1";
};

typedef std::vector< std::shared_ptr<COperatorResult> >           TOperatorResultList;
typedef std::map< std::string, std::shared_ptr< const COperatorResult> > TOperatorResultMap;

}

#endif

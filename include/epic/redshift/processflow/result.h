#ifndef _REDSHIFT_OPERATOR_RESULT_
#define _REDSHIFT_OPERATOR_RESULT_

#include <epic/core/common/datatypes.h>
#include <epic/core/common/managedobject.h>
#include <epic/core/common/constref.h>

#include <vector>
#include <ostream>
#include <map>

namespace NSEpic
{

class CDataStore;

class COperatorResult : public CManagedObject
{

public:

    COperatorResult();
    virtual ~COperatorResult();

    virtual Void Save( const CDataStore& store, std::ostream& stream ) const = 0;
    virtual Void SaveLine( const CDataStore& store, std::ostream& stream ) const = 0;
    //virtual Void Load( std::istream& stream ) = 0;

protected:

};

typedef std::vector< CConstRef<COperatorResult> >           TOperatorResultList;
typedef std::map< std::string, CConstRef<COperatorResult> > TOperatorResultMap;

}

#endif

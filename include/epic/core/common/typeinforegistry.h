#ifndef _CORE_COMMON_FACTORY_
#define _CORE_COMMON_FACTORY_

#include <vector>
#include <string>
#include <map>

#include <epic/core/common/datatypes.h>
#include <epic/core/common/singleton.h>
#include <epic/core/common/typeinfo.h>

#define TypeInfoRegistry (CTypeInfoRegistry::GetInstance())

namespace NSEpic
{

class CSerializable;

/**
 * \ingroup Core
 * Global registry of all type info for each managed object.
 *
 *
 */
class CTypeInfoRegistry : public CSingleton<CTypeInfoRegistry>
{

public:

    CTypeInfoRegistry();
    ~CTypeInfoRegistry();

    static Void             RegisterTypeInfo( CTypeInfo& objType );
    const CTypeInfo*        GetTypeInfo( String typeName ) const;

private:

    typedef std::map< std::string, CTypeInfo* > TypeInfoMap;
    TypeInfoMap                             m_TypeInfoMap;

    static CTypeInfo*                       m_RootEntry;

};

}

#endif

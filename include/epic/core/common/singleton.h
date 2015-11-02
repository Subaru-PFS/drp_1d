#ifndef _CORE_COMMON_SINGLETON_
#define _CORE_COMMON_SINGLETON_

#include <epic/core/common/datatypes.h>

namespace NSEpic
{

/**
 * \ingroup Core
 * Singleton base class
 */
template <typename T> 
class CSingleton
{

public:

    CSingleton( )
    {
        m_Instance = (T*) this;
    }

    virtual ~CSingleton( )
    {  
        m_Instance = NULL;  
    }

    static Bool IsCreated()
    {
        return m_Instance != NULL;
    }

    static T& GetInstance( )
    {  
        return *m_Instance ;
    }

private:

    static T*  m_Instance;
};

template <typename T> 
T* CSingleton<T>::m_Instance = NULL;


}

#endif

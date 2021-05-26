#ifndef _REDSHIFT_SPECTRUM_LSFFACTORY_
#define _REDSHIFT_SPECTRUM_LSFFACTORY_
/*
#include <map>
#include <stdio.h>
#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/spectrum/LSF.h"
#include "RedshiftLibrary/spectrum/LSF_NISPSIM_2016.h"
#include "RedshiftLibrary/spectrum/LSF_NISPVSSPSF_201707.h"
#include "RedshiftLibrary/spectrum/LSFConstantResolution.h"
#include "RedshiftLibrary/spectrum/LSFConstantWidth.h"
namespace NSEpic
{
class CLSF;

class CLSFFactory
{

public:

    template<typename...T>
    using CreateLSFFn = std::shared_ptr<CLSF> (*)(T&&...);//new way C++14 for typedef
    //static _CreateLSFFn CreateLSFFn; //define it as static to be called with no object instanciation
    
    ~CLSFFactory();

    static CLSFFactory *Get()
    {
        static CLSFFactory instance;
        return &instance;
    }

    void Register(const std::string &name, CreateLSFFn fn_makeLSF);
    
    template<typename... T> 
    std::shared_ptr<CLSF> Create(const std::string &name, T&&... args);//&& forward referencing

private:

    CLSFFactory();


    CLSFFactory(const CLSFFactory & other) = default; 
    CLSFFactory(CLSFFactory && other) = default; 
    CLSFFactory& operator=(const CLSFFactory& other) = default;  
    CLSFFactory& operator=(CLSFFactory&& other) = default; 

    typedef std::map<std::string, CreateLSFFn> FactoryMap;
    FactoryMap m_FactoryMap;
};
inline
CLSFFactory::CLSFFactory()
{
  Register("GaussianConstantWidth",      &CLSFGaussianConstantWidth::make_LSF);
  Register("GaussianConstantResolution", &CLSFGaussianConstantResolution::make_LSF);
  Register("GaussianNISPSIM2016",     &CLSFGaussianNISPSIM2016::make_LSF);
  Register("GaussianNISPVSSPSF201707",&CLSFGaussianNISPVSSPSF201707::make_LSF);
}

inline
CLSFFactory::~CLSFFactory()
{
  m_FactoryMap.clear();
}

inline
void CLSFFactory::Register(const std::string &name, CreateLSFFn fn_makeLSF)
{
    m_FactoryMap[name] = fn_makeLSF;
}

//there is a problem here cause we should pass variable nb of variables depending on the the selected object
template<typename... T>
std::shared_ptr<CLSF> CLSFFactory::Create(const std::string &name, T&&... args)
{
  FactoryMap::iterator it = m_FactoryMap.find(name);
  if(it!=m_FactoryMap.end())
    return it->second(std::forward<T>(args)...); 
  return NULL;
}

}*/
#endif

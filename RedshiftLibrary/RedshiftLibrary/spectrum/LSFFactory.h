#ifndef _REDSHIFT_SPECTRUM_LSFFACTORY_
#define _REDSHIFT_SPECTRUM_LSFFACTORY_
#define LSFFactory (CLSFFactory::GetInstance())

#include "RedshiftLibrary/common/singleton.h"

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
//the idea  is to define one single factory instance accross app instances
//now CLSFFactory, as defined here, "multiplexes" between the different LSFs
class CLSFFactory : public CSingleton<CLSFFactory>
{

public:

    using CreateLSFFn = std::shared_ptr<CLSF> (*)(const TLSFArguments& args);//new way C++14 for typedef

    static CLSFFactory *Get()
    {
        static CLSFFactory instance;
        return &instance;
    }

    void Register(const std::string &name, CreateLSFFn fn_makeLSF);
    
    std::shared_ptr<CLSF> Create(const std::string &name, TLSFArguments& args);

private:
    friend class CSingleton<CLSFFactory>;

    CLSFFactory();
    ~CLSFFactory();

    CLSFFactory(const CLSFFactory & other) = delete; 
    CLSFFactory(CLSFFactory && other) = delete; 
    CLSFFactory& operator=(const CLSFFactory& other) = delete;  
    CLSFFactory& operator=(CLSFFactory&& other) = delete; 

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

inline
std::shared_ptr<CLSF> CLSFFactory::Create(const std::string &name, TLSFArguments& args)
{
  FactoryMap::iterator it = m_FactoryMap.find(name);
  if(it!=m_FactoryMap.end())
    return it->second(args); 
  return NULL;
}

}
#endif

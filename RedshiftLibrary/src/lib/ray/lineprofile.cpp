#include "RedshiftLibrary/ray/lineprofile.h"
using namespace NSEpic;
using namespace std;

const std::string& CLineProfile::GetName(){
    return m_name;
}

//no need to define a constructor here

CLineProfile::CLineProfile(const Float64 nsigmasupport): 
m_nsigmasupport(nsigmasupport)
{}


CLineProfile::CLineProfile(const Float64 nsigmasupport, const std::string name): 
m_nsigmasupport(nsigmasupport),
m_name(name)
{}


Float64 CLineProfile::GetNSigmaSupport(){
    return m_nsigmasupport;
}
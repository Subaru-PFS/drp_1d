#ifndef _SOLVE_DESCRIPTION_H_
#define _SOLVE_DESCRIPTION_H_

#include <string>

namespace NSEpic
{


  class CSolveDescription
  {
  public:
    CSolveDescription(){}
    virtual ~CSolveDescription(){}
    
    static const std::string GetDescription(const std::string& method);
  };

}

#endif

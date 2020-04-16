#ifndef SETUP_FOR_FIND
#define SETUP_FOR_FIND
using namespace std;

namespace NSEpic
{


template< typename T >
class CFindExtrema
{

  public:
    
    CFindExtrema();
    ~CFindExtrema();
        
    void callFind( T chisquareResult );

};
    
#include <RedshiftLibrary/operator/findExtrema.hpp>
}

#endif  // SETUP_FOR_FIND

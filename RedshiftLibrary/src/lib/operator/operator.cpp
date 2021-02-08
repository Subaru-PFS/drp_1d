#include <RedshiftLibrary/operator/operator.h>

using namespace NSEpic;
using namespace std;

COperator::COperator()
{

}

COperator::~COperator()
{

}

Int32 COperator::getIndex(TFloat64List& list, Float64 z)
{
  std::vector<Float64>::iterator itr = std::lower_bound(list.begin(),list.end(), z);

  if (itr == list.end() || *itr != z)
    {
      size_t size = snprintf( nullptr, 0, "Could not find extrema solution index for %f", z) + 1; // Extra space for '\0'
      std::unique_ptr<char[]> buf( new char[ size ] );                                                        
      snprintf( buf.get(), size, "Could not find extrema solution index for %f", z);                          
      std::string _msg = std::string( buf.get(), buf.get() + size - 1 ); 
      throw runtime_error(_msg.c_str());
    }
  return (itr - list.begin()); 
}
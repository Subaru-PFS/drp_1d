#ifndef PYCONV_H
#define PYCONV_H
#include <RedshiftLibrary/common/datatypes.h>
namespace NSEpic
{
class PC
{
 public:
  static void get(const TFloat64List& vec,double ** data, int * size)
  {
    *data = const_cast<double*>(vec.data());
    *size = vec.size();
  }
  static void get(const TInt32List& vec,int ** data, int * size)
  {
    *data = const_cast<int*>(vec.data());
    *size = vec.size();
  }
};
}
#endif

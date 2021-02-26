#include <RedshiftLibrary/common/exception.h>

using namespace NSEpic;


AmzException::AmzException(ErrorCode ec,std::string message):
  _msg(message),
  code(ec)
{
  //std::ostringstream os;
  //os << boost::stacktrace::stacktrace();
  //stacktrace =  os.str();
  stacktrace = "";
}

AmzException::AmzException(const AmzException& e):
  _msg(e._msg),
  code(e.code),
  stacktrace(e.stacktrace)
{
}

const char* AmzException::getStackTrace() const {
  return stacktrace.c_str();
}


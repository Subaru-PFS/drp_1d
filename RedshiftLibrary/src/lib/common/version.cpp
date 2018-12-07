#include <RedshiftLibrary/version.h>

using namespace NSEpic;

const char* NSEpic::get_version() {
  static const char *version = CPF_REDSHIFT_REVISION;
  return version;
}

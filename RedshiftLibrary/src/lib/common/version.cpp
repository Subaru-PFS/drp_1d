#include <RedshiftLibrary/version.h>

using namespace NSEpic;

const char* NSEpic::get_version() {
  static char *version = CPF_REDSHIFT_REVISION;
  return version;
}

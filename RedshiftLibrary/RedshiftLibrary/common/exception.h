// ============================================================================
//
// This file is part of: AMAZED
//
// Copyright  Aix Marseille Univ, CNRS, CNES, LAM/CeSAM
//
// https://www.lam.fr/
//
// This software is a computer program whose purpose is to estimate the
// spectrocopic redshift of astronomical sources (galaxy/quasar/star)
// from there 1D spectrum.
//
// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use,
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info".
//
// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability.
//
// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or
// data to be ensured and,  more generally, to use and operate it in the
// same conditions as regards security.
//
// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.
// ============================================================================
#ifndef _REDSHIFT_EXCEPTION_
#define _REDSHIFT_EXCEPTION_

#include <exception>
#include <iostream>
#include <string>

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/log/log.h"

#define THROWG(code, msg)                                                      \
  throw AmzException(code, msg, __FILE__, __func__, __LINE__)
namespace NSEpic {

class AmzException : public std::exception {

public:
  AmzException(ErrorCode ec, const std::string &message, const char *filename_,
               const char *method_, int line_)
      : _msg(message), code(ec), filename(filename_), method(method_),
        line(line_) {}
  AmzException(const AmzException &e) = default;
  virtual ~AmzException() = default;
  AmzException(AmzException &&other) = default;
  AmzException &operator=(AmzException const &other) = default;
  AmzException &operator=(AmzException &&other) = default;

  virtual const char *what() const noexcept override { return _msg.c_str(); }

  const std::string &getMessage() const { return _msg; }

  ErrorCode getErrorCode() const { return code; }
  void setErrorCode(ErrorCode ec) { code = ec; }

  const std::string &getFileName() const { return filename; }
  void setFilename(std::string const &filename_) { filename = filename_; };
  const std::string &getMethod() const { return method; }
  void setMethod(std::string const &method_) { method = method_; };
  int getLine() const { return line; }
  void setLine(int line_) { line = line_; };
  void LogError(const std::string &msg = std::string()) const {
    Log.LogError(msg.empty() ? _msg : msg);
  };

protected:
  std::string _msg;
  ErrorCode code;
  std::string filename = "";
  std::string method = "";
  int line = -1;
};

} // namespace NSEpic

#endif

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
#ifndef _REDSHIFT_LOG_LOG_
#define _REDSHIFT_LOG_LOG_

#include "RedshiftLibrary/common/datatypes.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/mutex.h"
#include "RedshiftLibrary/common/singleton.h"

#include <cstdarg>

#define LOG_HANDLER_TABLE_SIZE 8
#define LOG_HANDLER_HEADER_LENGTH 64

#define Log (CLog::GetInstance())

namespace Log_test {
class consoleHandler_test;
class fileHandler_test;
} // namespace Log_test

namespace NSEpic {

class CLogHandler;
class AmzException;

/**
 * \ingroup Redshift
 * Responsible for logging features.
 * This class allows prioritized logging, and selective output of these log
 * messages according to handler's configuration.
 */
class CLog : public CSingleton<CLog> {

public:
  enum ELevel {
    nLevel_Critical = 100,
    nLevel_Error = 90,
    nLevel_Warning = 80,
    nLevel_Info = 70,
    nLevel_Detail = 65,
    nLevel_Debug = 60,
    nLevel_None = 0
  };

  friend class AmzException; // for calling private LogError

  // friend functions for calling private LogWarning :
  friend void CFlagWarning::warning(CFlagWarning::WarningCode c,
                                    std::string message);
  friend void CFlagWarning::warning(CFlagWarning::WarningCode c,
                                    const char *format, ...);

  void LogInfo(const char *format, ...);
  void LogDetail(const char *format, ...);
  void LogDebug(const char *format, ...);

  void LogInfo(const std::string &s);
  void LogDetail(const std::string &s);
  void LogDebug(const std::string &s);

  void log(const std::string &s, CLog::ELevel l);
  CMutex &GetSynchMutex();

  void Indent();
  void UnIndent();

private:
  friend class CSingleton<CLog>;
  friend class CLogHandler;

  friend class Log_test::consoleHandler_test;
  friend class Log_test::fileHandler_test;

  CLog();
  ~CLog();

  // LogError and LogWarning are private to be called from friend classes only
  // (CFlagWarning and AmzException)
  void LogError(const char *format, ...);
  void LogWarning(const char *format, ...);

  void LogError(const std::string &s);
  void LogWarning(const std::string &s);

  void LogEntry(ELevel lvl, const char *format, va_list &args);
  void AddHandler(CLogHandler &handler);
  void RemoveHandler(CLogHandler &handler);
  const char *GetHeader(CLog::ELevel lvl);

  // Attributes
  CLogHandler *m_HandlerTable[LOG_HANDLER_TABLE_SIZE];
  Char m_CurrentHeader[LOG_HANDLER_HEADER_LENGTH];
  Char *m_WorkingBuffer;
  Char m_IndentBuffer[128];
  Int32 m_IndentCount;

  CMutex m_Mutex;
};

} // namespace NSEpic

#endif

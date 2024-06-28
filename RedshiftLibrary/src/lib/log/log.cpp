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
#include <chrono>

#include "RedshiftLibrary/common/defaults.h"
#include "RedshiftLibrary/log/handler.h"
#include "RedshiftLibrary/log/log.h"

using namespace NSEpic;

/**
 * Creates a singleton that holds a buffer and a table of the log handlers.
 */
CLog::CLog() {
  Int32 i;
  for (i = 0; i < LOG_HANDLER_TABLE_SIZE; i++) {
    m_HandlerTable[i] = NULL;
  }
}

CLog::~CLog() {}

void CLog::logEntry(const std::string &msg, CLog::ELevel logLevel,
                    bool withTimestamp) {
  m_Mutex.Lock();
  std::string messageToLog = GetHeader(logLevel) + msg;
  if (withTimestamp)
    messageToLog = timeStampString() + messageToLog;
  for (int i = 0; i < LOG_HANDLER_TABLE_SIZE; i++) {
    if (m_HandlerTable[i]) {
      if (m_HandlerTable[i]->GetLevelMask() <= logLevel) {
        m_HandlerTable[i]->LogEntry(logLevel, messageToLog);
      }
    }
  }
  m_Mutex.Unlock();
}

std::string CLog::timeStampString() {
  // Get the current system time point
  auto now = std::chrono::system_clock::now();

  // Convert the time point to a time_t object (C-style time)
  std::time_t time = std::chrono::system_clock::to_time_t(now);

  // Convert to string
  std::string timeString = std::ctime(&time);

  // Remove last character (newline)
  timeString = timeString.erase(timeString.size() - 1) + " ";

  return timeString;
}

void CLog::LogError(const std::string &msg, bool withTimestamp) {
  logEntry(msg, nLevel_Error, withTimestamp);
}
void CLog::LogDebug(const std::string &msg, bool withTimestamp) {
  logEntry(msg, nLevel_Debug, withTimestamp);
}
void CLog::LogInfo(const std::string &msg, bool withTimestamp) {
  logEntry(msg, nLevel_Info, withTimestamp);
}
void CLog::LogWarning(const std::string &msg, bool withTimestamp) {
  logEntry(msg, nLevel_Warning, withTimestamp);
}
void CLog::LogDetail(const std::string &msg, bool withTimestamp) {
  logEntry(msg, nLevel_Detail, withTimestamp);
}

/**
 * Remove the reference for the input handler if it exists in the handler table.
 */
void CLog::RemoveHandler(CLogHandler &handler) {
  Int32 i;

  for (i = 0; i < LOG_HANDLER_TABLE_SIZE; i++) {
    if (m_HandlerTable[i] == &handler) {
      m_HandlerTable[i] = NULL;
      break;
    }
  }
}

/**
 * Insert a reference to the input handler in the first vacant position in the
 * handler table. Will _not_ prevent a handler to be included multiple times.
 */
void CLog::AddHandler(CLogHandler &handler) {
  Int32 i;

  for (i = 0; i < LOG_HANDLER_TABLE_SIZE; i++) {
    if (m_HandlerTable[i] == NULL) {
      m_HandlerTable[i] = &handler;
      break;
    }
  }
}

CMutex &CLog::GetSynchMutex() { return m_Mutex; }

/**
 * Given a message priority level, returns an appropriate prefix to the message.
 */
std::string CLog::GetHeader(CLog::ELevel logLevel) {
  std::string header;
  switch (logLevel) {
  case nLevel_Critical:
    header = "Critical: ";
    break;
  case nLevel_Error:
    header = "Error: ";
    break;
  case nLevel_Warning:
    header = "Warning: ";
    break;
  case nLevel_Info:
    header = "Info: ";
    break;
  case nLevel_Detail:
    header = "Detail: ";
    break;
  case nLevel_Debug:
    header = "Debug: ";
    break;
  case nLevel_None:
    header = "";
    break;
  default:
    return NULL;
  }

  return header;
}

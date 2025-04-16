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
#include <fstream>
#include <iostream>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <regex>

#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/log/consolehandler.h"
#include "RedshiftLibrary/log/filehandler.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Log_test)

namespace bfs = boost::filesystem;

std::pair<std::string, std::string> splitString(int n,
                                                const std::string &input) {
  if (ssize(input) <= n) {
    return {input, ""};
  } else {
    return {input.substr(0, n), input.substr(n)};
  }
}

BOOST_AUTO_TEST_CASE(consoleHandler_test) {

  // Add consoleHandler to the log
  CLogConsoleHandler console_handler;

  // check default level
  BOOST_CHECK(console_handler.GetLevelMask() == Log.nLevel_Warning);

  // Create one message
  std::string message(Formatter() << "test log " << 1);

  // Log
  Log.LogWarning(message);
  console_handler.LogEntry(80, "Warning" + message);

  console_handler.SetLevelMask(65);
  BOOST_CHECK(console_handler.GetLevelMask() == Log.nLevel_Detail);
}

BOOST_AUTO_TEST_CASE(fileHandler_test) {

  bfs::path logFile = bfs::unique_path("log_%%%%%%%%%%");

  // Add consoleHandler to the log
  CLogFileHandler file_handler(logFile.c_str());

  // check default level
  BOOST_CHECK(file_handler.GetLevelMask() == Log.nLevel_Warning);

  file_handler.SetLevelMask(0);
  BOOST_CHECK(file_handler.GetLevelMask() == Log.nLevel_None);

  std::string message;
  message = Formatter() << "test log warning";
  Log.LogWarning(message);

  message = Formatter() << "test log info";
  Log.LogInfo(message);

  message = Formatter() << "test log error";
  Log.LogError(message);

  message = Formatter() << "test log detail";
  Log.LogDetail(message);

  message = Formatter() << "test log debug";
  Log.LogDebug(message);

  message = Formatter() << "test log none";
  Log.logEntry(message, Log.nLevel_None);

  message = Formatter() << "test log critical";
  Log.logEntry(message, Log.nLevel_Critical);

  message = Formatter() << "test log warning";
  file_handler.LogEntry(80, "Warning: " + message);

  message = Formatter() << "test log warning without timestamp";
  Log.logEntry(message, Log.nLevel_Warning, false);

  CMutex &mutex = Log.GetSynchMutex();

  std::string line;
  std::ifstream myfile(logFile.c_str());
  TStringList lines;

  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      lines.push_back(line);
    }
    myfile.close();
  }

  std::string line_ref;
  std::string timestamp;
  std::string logMessage;
  std::regex timestampPattern(
      "\\b\\w{3} \\w{3}\\s{1,2}\\d{1,2} \\d{2}:\\d{2}:\\d{2} \\d{4}\\b");
  for (Int32 i = 0; i < ssize(lines); i++) {
    if (i == 0)
      line_ref = "Warning: test log warning";
    else if (i == 1)
      line_ref = "Info: test log info";
    else if (i == 2)
      line_ref = "Error: test log error";
    else if (i == 3)
      line_ref = "Detail: test log detail";
    else if (i == 4)
      line_ref = "Debug: test log debug";
    else if (i == 5)
      line_ref = "test log none";
    else if (i == 6)
      line_ref = "Critical: test log critical";
    else if (i == 7)
      line_ref = "Warning: test log warning";
    else if (i == 8)
      line_ref = "Warning: test log warning without timestamp";

    if (i < 7) {
      // Remove all timestamps for logs containing one
      auto parts = splitString(25, lines[i]);
      timestamp = parts.first;
      logMessage = parts.second;
      BOOST_CHECK(std::regex_search(timestamp, timestampPattern));
    } else {
      logMessage = lines[i];
    }
    BOOST_CHECK_EQUAL(logMessage, line_ref);
  }

  bfs::remove_all(logFile);
}

BOOST_AUTO_TEST_SUITE_END()
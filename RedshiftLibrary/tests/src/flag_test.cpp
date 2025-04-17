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
#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/common/formatter.h"
#include "RedshiftLibrary/common/size.h"
#include "RedshiftLibrary/log/consolehandler.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(Flag_test)

BOOST_AUTO_TEST_CASE(some_test) {

  TWarningMsgList listMsg_ref;

  // Add consoleHandler to the log
  CLogConsoleHandler console_handler;

  // Create one warning with digit 3
  std::string message(Formatter()
                      << "test log warning : code "
                      << (int)WarningCode::AIR_VACUUM_CONVERSION_IGNORED);
  listMsg_ref.push_back(
      std::make_pair(WarningCode::AIR_VACUUM_CONVERSION_IGNORED, message));
  Flag.warning(WarningCode::AIR_VACUUM_CONVERSION_IGNORED, message);

  // Recover flagList
  TWarningMsgList listMsg = Flag.getListMessages();
  for (Int32 i = 0; i < ssize(listMsg); i++) {
    BOOST_CHECK(listMsg[i].first == listMsg_ref[i].first);
    BOOST_CHECK(listMsg[i].second == listMsg_ref[i].second);
  }

  // Create one warning with digit 11
  message = Formatter()
            << "test log warning : code "
            << (int)WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM;
  listMsg_ref.push_back(std::make_pair(
      WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM, message));
  Flag.warning(WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM,
               "test log warning : code %d",
               WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM);

  // Recover flagList
  listMsg = Flag.getListMessages();
  for (Int32 i = 0; i < ssize(listMsg); i++) {
    BOOST_CHECK(listMsg[i].first == listMsg_ref[i].first);
    BOOST_CHECK(listMsg[i].second == listMsg_ref[i].second);
  }

  Int32 refVal =
      1 << int(WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM) |
      1 << int(WarningCode::AIR_VACUUM_CONVERSION_IGNORED);

  // Recover flag value
  BOOST_CHECK(Flag.getBitMask() == refVal); // 2056 or 00001000 00001000

  // Create one warning with another digit 12
  message = Formatter()
            << "test log warning 2 : code "
            << (int)WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM;
  listMsg_ref.push_back(std::make_pair(
      WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM, message));
  Flag.warning(WarningCode::FORCED_CONTINUUM_COMPONENT_TO_FROMSPECTRUM,
               message);

  // Recover flagList
  listMsg = Flag.getListMessages();
  for (Int32 i = 0; i < ssize(listMsg); i++) {
    BOOST_CHECK(listMsg[i].first == listMsg_ref[i].first);
    BOOST_CHECK(listMsg[i].second == listMsg_ref[i].second);
  }

  // Reset mask
  Flag.resetFlag();
  listMsg = Flag.getListMessages();
  BOOST_CHECK(Flag.getBitMask() == 0);
  BOOST_CHECK(listMsg.size() == 0);
}

BOOST_AUTO_TEST_SUITE_END()

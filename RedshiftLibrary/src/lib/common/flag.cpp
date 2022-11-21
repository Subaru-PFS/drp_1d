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
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/log/log.h"

#define FLAG_MSG_WORKING_BUFFER_SIZE 4096

using namespace NSEpic;

CFlagWarning::CFlagWarning() {
  m_WorkingBuffer = new Char[FLAG_MSG_WORKING_BUFFER_SIZE];
}

CFlagWarning::~CFlagWarning() { delete[] m_WorkingBuffer; }

void CFlagWarning::warning(WarningCode c, std::string message) {
  m_flags |= (1 << (int)c);
  m_messageList.push_back(make_pair(c, message));
  Log.LogWarning(message);
}

void CFlagWarning::warning(WarningCode c, const char *format, ...) {
  m_flags |= (1 << (int)c);

  va_list args;
  va_start(args, format);
  vsnprintf(m_WorkingBuffer, FLAG_MSG_WORKING_BUFFER_SIZE, format, args);
  va_end(args);

  std::string message(m_WorkingBuffer);
  m_messageList.push_back(make_pair(c, message));
  Log.LogWarning(message);
}

Int32 CFlagWarning::getBitMask() { return m_flags; }

void CFlagWarning::resetFlag() {
  m_flags = 0;
  m_messageList.clear();
}

const TWarningMsgList &CFlagWarning::getListMessages() { return m_messageList; }

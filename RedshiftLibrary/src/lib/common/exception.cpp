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
#include "RedshiftLibrary/common/exception.h"
#include "RedshiftLibrary/common/flag.h"
#include "RedshiftLibrary/operator/flagResult.h"
#include "RedshiftLibrary/processflow/context.h"
#include "RedshiftLibrary/processflow/resultstore.h"

using namespace NSEpic;

GlobalException::GlobalException(ErrorCode ec, const std::string &message,
                                 const char *filename_, const char *method_,
                                 int line_) noexcept
    : AmzException(ec, message, filename_, method_, line_) {

  std::shared_ptr<COperatorResultStore> resultStore = Context.GetResultStore();

  if (resultStore->getScopeDepth() >= 2) {
    if (!resultStore->hasCurrentMethodWarningFlag()) {
      resultStore->StoreGlobalResult(
          resultStore->GetScopedNameAt("warningFlag", 2),
          std::make_shared<const CFlagLogResult>(Flag.getBitMask(),
                                                 Flag.getListMessages()));
    } else
      LogError(Formatter() << "Warning flag already exists for "
                           << resultStore->getCurrentScopeNameAt(2));
  } else {
    if (!resultStore->hasContextWarningFlag()) {
      resultStore->StoreGlobalResult(
          "context_warningFlag",
          std::make_shared<const CFlagLogResult>(Flag.getBitMask(),
                                                 Flag.getListMessages()));
    } else
      LogError("Context warning flag already exists");
  }
  Flag.resetFlag();
}

InternalException::InternalException(ErrorCode ec, const std::string &message,
                                     const char *filename_, const char *method_,
                                     int line_) noexcept
    : AmzException(ec, message, filename_, method_, line_) {}

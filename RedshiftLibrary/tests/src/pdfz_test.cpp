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
#include <regex>

#include <boost/test/unit_test.hpp>

#include "RedshiftLibrary/operator/pdfz.h"

using namespace NSEpic;
using namespace std;

BOOST_AUTO_TEST_SUITE(Pdfz_test)

BOOST_AUTO_TEST_CASE(checkWindowSize_test) {
  COperatorPdfz op_pdfz("something",
                        0.0, // no peak Separation in 2nd pass
                        0.0, // cut threshold
                        5,   // max nb of final (2nd pass) candidates
                        true,
                        "SPE", // Id_prefix
                        false, // do not allow extrema at border
                        1      // one peak/window only
  );

  TFloat64Range integration_range;
  TFloat64Range window_range;

  integration_range = TFloat64Range{0.1, 0.2};
  window_range = TFloat64Range{0.05, 0.15};
  op_pdfz.checkWindowSize(integration_range, window_range);
  BOOST_CHECK(Flag.getListMessages().size() == 1);
  BOOST_CHECK(Flag.getListMessages()[0].first ==
              WarningCode::PDF_INTEGRATION_WINDOW_TOO_SMALL);

  Flag.resetFlag();
  integration_range = TFloat64Range{0.05, 0.15};
  window_range = TFloat64Range{0.1, 0.2};
  op_pdfz.checkWindowSize(integration_range, window_range);
  BOOST_CHECK(Flag.getListMessages().size() == 1);
  BOOST_CHECK(Flag.getListMessages()[0].first ==
              WarningCode::PDF_INTEGRATION_WINDOW_TOO_SMALL);

  Flag.resetFlag();
  integration_range = TFloat64Range{0.1, 0.2};
  window_range = TFloat64Range{0.05, 0.25};
  op_pdfz.checkWindowSize(integration_range, window_range);
  BOOST_CHECK(Flag.getListMessages().size() == 0);
}

BOOST_AUTO_TEST_SUITE_END()

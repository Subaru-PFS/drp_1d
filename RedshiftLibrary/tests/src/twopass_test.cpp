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
#include "RedshiftLibrary/common/indexing.h"
#include "RedshiftLibrary/operator/twopass.h"
#include <boost/test/unit_test.hpp>
using namespace NSEpic;
BOOST_AUTO_TEST_SUITE(TwoPass)

BOOST_AUTO_TEST_CASE(spanRedshift_test) {
  Float64 z = 5.;
  Float64 step = log(1 + 0.5);
  TFloat64List redshifts{0, 5, 9};
  bool zLogSampling = true;
  Float64 secondPass_halfwindowsize = log(1 + 1);
  Int32 ref_idx = 1;
  TFloat64List extendedRedshifts_ref{3, 5, 8};

  // prepare object
  COperatorTwoPass op;
  op.m_Redshifts = redshifts;
  op.m_fineStep = step;
  op.m_zLogSampling = zLogSampling;
  op.m_secondPass_halfwindowsize = secondPass_halfwindowsize;

  TFloat64List extendedList = op.spanRedshiftWindow(z);
  // check is sorted
  BOOST_CHECK(
      std::is_sorted(std::begin(extendedList), std::end(extendedList)) == true);
  BOOST_CHECK(extendedList == extendedRedshifts_ref);
  // check presence of z in extendedList
  Int32 idx = CIndexing<Float64>::getIndex(extendedList, z);
  BOOST_CHECK(idx == ref_idx);
}

BOOST_AUTO_TEST_SUITE_END()

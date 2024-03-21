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

#include "RedshiftLibrary/processflow/result.h"

using namespace NSEpic;

BOOST_AUTO_TEST_SUITE(operatorResult)

BOOST_AUTO_TEST_CASE(operatorResult_test) {
  COperatorResult result("COperatorResult");
  std::string type = result.getType();

  BOOST_CHECK(type == "COperatorResult");

  // virtual function
  BOOST_CHECK_THROW(result.getCandidateDatasetType("dataset"), AmzException);
  BOOST_CHECK_THROW(result.getCandidate(0, "dataset"), AmzException);
  BOOST_CHECK_THROW(result.HasCandidateDataset("dataset"), AmzException);

  // copy and copy assignement
  COperatorResult result2(result);
  type = result2.getType();
  BOOST_CHECK(type == "COperatorResult");

  COperatorResult result3("COperatorResult3");
  result3 = result2;
  type = result3.getType();
  BOOST_CHECK(type == "COperatorResult");

  // move and move assignement
  COperatorResult result4(std::move(result));
  type = result4.getType();
  BOOST_CHECK(type == "COperatorResult");

  COperatorResult result5("COperatorResult5");
  result5 = std::move(result2);
  type = result5.getType();
  BOOST_CHECK(type == "COperatorResult");
}

BOOST_AUTO_TEST_SUITE_END()
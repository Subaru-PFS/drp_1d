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
#include "RedshiftLibrary/line/catalog.h"
#include "RedshiftLibrary/line/line.h"
#include "RedshiftLibrary/operator/linematching.h"

#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace NSEpic;
BOOST_AUTO_TEST_SUITE(test_linematching)

Float64 precision = 1e-12;
BOOST_AUTO_TEST_CASE(Compute) {

  CLineDetectedCatalog detectedLineCatalog;
  CLineProfile_ptr profilesym{
      std::unique_ptr<CLineProfileSYM>(new CLineProfileSYM())};
  detectedLineCatalog.Add(CLineDetected(
      "Line1", 20, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));
  detectedLineCatalog.Add(CLineDetected(
      "Line2", 40, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));

  detectedLineCatalog.Add(CLineDetected(
      "Line3", 250, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));
  detectedLineCatalog.Add(CLineDetected(
      "Line4", 800, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));

  CLineCatalog restLineCatalog;
  restLineCatalog.Add(CLine("Line2_1", 10, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 0, "0"));
  restLineCatalog.Add(CLine("Line2_2", 20, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 1, "1"));
  restLineCatalog.Add(CLine("Line2_3", 125, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 2, "2"));
  restLineCatalog.Add(CLine("Line2_4", 400, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 3, "3"));

  CLineMatching lineMathing = CLineMatching();
  TFloat64Range redshiftRange = TFloat64Range(0, 5);
  Int32 nThreshold = 4;
  Float64 tol = 0.05;
  auto typeFilter = CLine::EType::nType_Emission;
  auto detectedForceFilter = CLine::EForce::nForce_Strong;
  auto restRorceFilter = CLine::EForce::nForce_Strong;

  std::shared_ptr<CLineMatchingResult> res = lineMathing.Compute(
      detectedLineCatalog, restLineCatalog, redshiftRange, nThreshold, tol,
      typeFilter, detectedForceFilter, restRorceFilter);
  BOOST_CHECK_EQUAL(res->SolutionSetList.size(), 1);

  // BOOST_CHECK_EQUAL(res->SolutionSetList[0].size(), 1);
  BOOST_CHECK_CLOSE(res->SolutionSetList[0][0].Redshift, 1, precision);

  detectedLineCatalog = CLineDetectedCatalog();
  detectedLineCatalog.Add(CLineDetected(
      "Line1", 800, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));

  restLineCatalog = CLineCatalog();
  restLineCatalog.Add(CLine("Line2_1", 10, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 0, "0"));
  restLineCatalog.Add(CLine("Line2_2", 20, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 1, "1"));
  restLineCatalog.Add(CLine("Line2_3", 125, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 2, "2"));
  restLineCatalog.Add(CLine("Line2_4", 400, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 3, "3"));
  nThreshold = 1;
  res = lineMathing.Compute(detectedLineCatalog, restLineCatalog, redshiftRange,
                            nThreshold, tol, typeFilter, detectedForceFilter,
                            restRorceFilter);
  BOOST_CHECK_EQUAL(res->SolutionSetList.size(), 1);

  // BOOST_CHECK_EQUAL(res->SolutionSetList[0].size(), 1);
  BOOST_CHECK_CLOSE(res->SolutionSetList[0][0].Redshift, 1, precision);

  detectedLineCatalog = CLineDetectedCatalog();
  detectedLineCatalog.Add(CLineDetected(
      "Line1", 20, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));
  detectedLineCatalog.Add(CLineDetected(
      "Line2", 20, CLine::EType::nType_Absorption, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));
  detectedLineCatalog.Add(CLineDetected(
      "Line3", 40, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));
  detectedLineCatalog.Add(CLineDetected(
      "Line4", 80, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));
  detectedLineCatalog.Add(CLineDetected(
      "Line5", 160, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));
  detectedLineCatalog.Add(CLineDetected(
      "Line6", 320, CLine::EType::nType_Emission, profilesym->Clone(),
      CLine::EForce::nForce_Strong, 1, 4, 5.6));

  restLineCatalog = CLineCatalog();
  restLineCatalog.Add(CLine("Line2_1", 10, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 0, "0"));
  restLineCatalog.Add(CLine("Line2_2", 10, CLine::EType::nType_Absorption,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 1, "1"));
  restLineCatalog.Add(CLine("Line2_3", 20, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 2, "2"));
  restLineCatalog.Add(CLine("Line2_4", 40, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 3, "3"));
  restLineCatalog.Add(CLine("Line2_5", 80, CLine::EType::nType_Emission,
                            profilesym->Clone(), CLine::EForce::nForce_Strong,
                            0., false, "-1", 1.0, "em", 4, "4"));
  nThreshold = 3;
  redshiftRange = TFloat64Range(0, 10);
  /*  res = lineMathing.Compute(detectedLineCatalog, restLineCatalog,
  redshiftRange, nThreshold, tol, typeFilter, detectedForceFilter,
  restRorceFilter); BOOST_CHECK_EQUAL(res->SolutionSetList.size(), 3);

  //BOOST_CHECK_EQUAL(res->SolutionSetList[0].size(), 1);
  BOOST_CHECK_CLOSE(res->SolutionSetList[0][0].Redshift,1, precision);
  BOOST_CHECK_CLOSE(res->SolutionSetList[1][0].Redshift,3, precision);
  BOOST_CHECK_CLOSE(res->SolutionSetList[2][0].Redshift,7, precision);
  */
}

BOOST_AUTO_TEST_SUITE_END()
